#include "sound.hpp"
#include <cstring>
#include <algorithm>
#include <cmath>
#include <stdexcept>

#define CLIP(x,min,max) ((x)<(min)?(min):((x)>(max)?(max):(x)))

S3W2_Sound::S3W2_Sound()
{
    initialize();
}

S3W2_Sound::~S3W2_Sound()
{
    pcm_ram.reset();
}

void S3W2_Sound::initialize()
{
    // allocate PCM RAM
    pcm_ram = std::make_unique<std::array<uint8_t, PCM_RAM_SIZE>>();
    std::fill(pcm_ram->begin(), pcm_ram->end(), 0);
    for (int i = 0; i < NUM_CHANNELS; ++i)
        resetChannel(i);
}

void S3W2_Sound::reset()
{
    if (pcm_ram)
        std::fill(pcm_ram->begin(), pcm_ram->end(), 0);
    for (int i = 0; i < NUM_CHANNELS; ++i)
        resetChannel(i);
}

void S3W2_Sound::resetChannel(int ch)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return;
    Channel &c = channels[ch];
    c.frequency = 0;
    c.waveform_type = WAVEFORM_WAVETABLE;
    c.volume = 0;
    c.panpot = 0xFF;
    c.modulation_type = MOD_NONE;
    c.modulation_param_1 = 0;
    c.modulation_param_2 = 0;
    c.modulation_target = 0;
    c.modulation_targeting_mode = 0;
    c.pcm_start_addr = 0;
    c.pcm_end_addr = 0;
    c.pcm_loop_addr = 0;
    c.pcm_control = 0;
    c.sample_control = 0;
    c.phase = 0.0;
    c.sqv4_ready = false;
    c.old_phase = 0;
    c.slew_cur = 0;
    c.last_sample = 0;
    c.lfsr_state = 0x12D4803C;
    c.active = true;
    c.wavetable.fill(0x80);
    resetFIFO(c);
}

void S3W2_Sound::resetFIFO(Channel& c) {
    c.fifo.clear();
    c.fifo_held=0;
    c.fifo_hold=0;
    c.fifo_sub=0;
    c.fifo_latch=0;
}

void S3W2_Sound::setRevision(int revision) {
    bool nextA=revision>=1;
    bool nextB=revision>=2;
    if (nextA==revisionA && nextB==revisionB) return;
    revisionA=nextA;
    revisionB=nextB;
    for (auto& c: channels) {
        resetFIFO(c);
        if (!revisionB) c.sample_control=0;
        c.sqv4_ready=false;
        if (!revisionA && c.waveform_type==WAVEFORM_FIFO_DAC)
            c.waveform_type=WAVEFORM_WAVETABLE;
    }
}

void S3W2_Sound::writeRegister(uint32_t address, uint8_t value)
{
    // basic routing: wavetable area 0x000 - 0x7FF, control area 0x800 - 0x8FF
    if (address <= 0x7FF)
    {
        // each channel has 256 bytes, 8 channels => 0x000-0x7FF
        uint32_t ch = (address >> 8) & 0x7; // address / 256
        uint8_t offset = address & 0xFF;
        writeChannelWavetable(static_cast<int>(ch), offset, value);
        return;
    }
    else if (address >= 0x800 && address <= 0x8FF)
    {
        uint32_t ch = ((address - 0x800) >> 5) & 0x7; // each channel 32 bytes
        uint8_t offset = (address - 0x800) & 0x1F;
        if (offset < 0x10) {
            writeChannelControl(static_cast<int>(ch), offset, value);
        } else {
            writeChannelPCMReg(static_cast<int>(ch), offset - 0x10, value);
        }
        return;
    }
    // other addresses not implemented yet
    (void)address;
    (void)value;
}

uint8_t S3W2_Sound::readRegister(uint32_t address)
{
    if (address <= 0x7FF)
    {
        uint32_t ch = (address >> 8) & 0x7;
        uint8_t offset = address & 0xFF;
        return readChannelWavetable(static_cast<int>(ch), offset);
    }
    else if (address >= 0x800 && address <= 0x8FF)
    {
        uint32_t ch = ((address - 0x800) >> 5) & 0x7;
        uint8_t offset = (address - 0x800) & 0x1F;
        return readChannelControl(static_cast<int>(ch), offset);
    }
    return 0;
}

std::vector<std::vector<std::vector<int16_t>>> S3W2_Sound::clock(size_t samples)
{
  // Prepare output buffer: [channel][left/right][sample]
    std::vector<std::vector<std::vector<int16_t>>> output(NUM_CHANNELS,
        std::vector<std::vector<int16_t>>(2, std::vector<int16_t>(samples, 0)));
    for (size_t i = 0; i < samples; ++i)
    {
        for (int ch = 0; ch < NUM_CHANNELS; ++ch)
        {
            int32_t left=0, right=0;
            for (int j = 0; j < 4; ++j)
            {
                int16_t sample = generateSample(ch);
                uint8_t panL = channels[ch].panpot >> 4; // upper 4 bits for left
                uint8_t panR = channels[ch].panpot & 0x0F; // lower 4 bits for right
                left += sample * panL / 15;
                right += sample * panR / 15;
                //printf("Ch %d Sample %d, Subsample %d: %d\n", ch, i, j, sample);
            }
            // FIFO samples already have full 16-bit gain. Average the four
            // internal ticks without overflowing the channel output buffer.
            bool fifo=revisionA && channels[ch].waveform_type==WAVEFORM_FIFO_DAC;
            output[ch][0][i]=static_cast<int16_t>(fifo?left/4:left);
            output[ch][1][i]=static_cast<int16_t>(fifo?right/4:right);
        }
    }
    return output;
}

int16_t S3W2_Sound::generateSample(int ch)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    Channel &c = channels[ch];
    if (!c.active)
        return 0;
    if (revisionA && c.waveform_type==WAVEFORM_FIFO_DAC) {
        int16_t sample=generateFIFOSample(ch);
        c.last_sample=sample>>8;
        return sample; // fixed gain; VOLUME does not apply
    }
    if (revisionA && (c.waveform_type==WAVEFORM_WAVETABLE || c.waveform_type==WAVEFORM_PCM)) {
        int32_t q8=generateInterpolatedSample(ch,c.waveform_type==WAVEFORM_PCM);
        c.last_sample=static_cast<int16_t>(q8>>8);
        return static_cast<int16_t>(q8*c.volume/(4*256));
    }
    if (c.waveform_type==WAVEFORM_SQV4) {
        int32_t q8=generateSQV4Sample(ch);
        c.last_sample=static_cast<int16_t>(q8>>8);
        return static_cast<int16_t>(q8*c.volume/(255*256*4));
    }
    int16_t sample = 0;
    switch (c.waveform_type)
    {
    case WAVEFORM_WAVETABLE:
        sample = generateWavetableSample(ch);
        break;
    case WAVEFORM_PCM:
        sample = generatePCMSample(ch);
        break;
    case WAVEFORM_NOISE:
        sample = generateNoiseSample(ch);
        break;
    case WAVEFORM_DMA_PCM:
        sample = generateDMAPCMSample(ch);
        break;
    default:
        sample = 0;
        break;
    }
    c.last_sample = sample;
    return sample*c.volume/4;
}

int16_t S3W2_Sound::generateFIFOSample(int ch) {
    Channel& c=channels[ch];
    if (c.fifo_sub==0) {
        if (c.fifo_hold) {
            --c.fifo_hold;
        } else if (!c.fifo.empty()) {
            c.fifo_held=c.fifo.front();
            c.fifo.pop_front();
            c.fifo_hold=c.frequency;
        }
    }
    c.fifo_sub=(c.fifo_sub+1)&3;
    return c.fifo_held;
}

// Unsigned Q8.8 interpolation preserves fractional steps through modulation
// and volume. Modulation IDs retain the Furnace core's existing mapping.
int32_t S3W2_Sound::generateInterpolatedSample(int ch, bool pcm) {
    Channel& c=channels[ch];
    const uint32_t waveMask=(256u>>(c.sample_control&3))-1;
    // The chip compensates for the table length: the same FREQUENCY gives the
    // same pitch for 256/128/64/32-sample tables.
    if (!pcm) c.phase+=static_cast<double>(c.frequency)/SOUND_CLOCK*(waveMask+1)/256.0;
    double scaled=c.phase*256;
    uint64_t phase=static_cast<uint64_t>(scaled);
    int frac=CLIP(static_cast<int>((scaled-phase)*256),0,255);
    if (c.sample_control&4) frac=0; // interpolation off: hold the integer position
    int target=convertToAbsoluteChannelAddress(ch,c.modulation_targeting_mode,c.modulation_target);
    if (c.modulation_type==MOD_PHASE) {
        int64_t offset=(c.modulation_param_1*channels[target].last_sample)>>12;
        phase+=offset;
    }
    uint32_t addr=pcm?((phase+c.pcm_start_addr)&(PCM_RAM_SIZE-1)):(phase&waveMask);
    uint32_t next;
    int current, following;
    if (pcm) {
        next=addr+1;
        if (next>=c.pcm_end_addr)
            next=(c.pcm_control&2)?c.pcm_loop_addr:addr;
        current=(*pcm_ram)[addr];
        following=(*pcm_ram)[next&(PCM_RAM_SIZE-1)];
    } else {
        current=c.wavetable[addr];
        following=c.wavetable[(addr+1)&waveMask];
    }
    int32_t q8=current*(256-frac)+following*frac;
    if (c.modulation_type==MOD_RING) {
        int64_t ring=static_cast<int64_t>(q8-(128<<8))*channels[target].last_sample*c.modulation_param_1/65536+(128<<8);
        q8=static_cast<int32_t>(CLIP(ring,0,255<<8));
    } else if (c.modulation_type==MOD_SLEW_RATE || c.modulation_type==MOD_SLEW_STEP) {
        int diff=q8-c.slew_cur;
        int step=c.modulation_param_1;
        c.slew_cur+=c.modulation_type==MOD_SLEW_RATE?CLIP(diff,-step,step):(diff>0?step:-step);
        c.slew_cur=CLIP(c.slew_cur,0,255<<8);
        q8=c.slew_cur;
    }
    int32_t out=q8-(128<<8);
    if (pcm) {
        if (addr>=c.pcm_end_addr) {
            if (c.pcm_control&2) {
                c.phase=std::max(0,static_cast<int>(c.pcm_loop_addr)-static_cast<int>(c.pcm_start_addr))/256.0;
            } else {
                c.phase=std::max(0,static_cast<int>(c.pcm_end_addr)-static_cast<int>(c.pcm_start_addr))/256.0;
                out=0;
            }
        } else {
            c.phase+=static_cast<double>(c.frequency)/SOUND_CLOCK/8;
        }
    }
    return out;
}

int16_t S3W2_Sound::generateWavetableSample(int ch)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    Channel &c = channels[ch];
    if (c.waveform_type != WAVEFORM_WAVETABLE)
        return 0;
    // calculate phase increment
    c.phase+=static_cast<double>(c.frequency)/SOUND_CLOCK; // advance phase
    uint8_t index = 0;
    int16_t phase_offset = 0;
    uint64_t phase = static_cast<uint64_t>(c.phase*256) & 0xFFFFFFFFFFFFFFFF;
    uint8_t abs_target_ch = convertToAbsoluteChannelAddress(static_cast<uint8_t>(ch), c.modulation_targeting_mode, c.modulation_target);
    switch (c.modulation_type)
    {
    case MOD_NONE:
        index = static_cast<uint8_t>(phase & 0xFF);
        break;
    case MOD_PHASE:
        phase_offset = static_cast<int16_t>((c.modulation_param_1 * channels[abs_target_ch].last_sample) >> 12);
        index = static_cast<uint8_t>((phase + phase_offset) & 0xFF);
        break;
    default:
        index = static_cast<uint8_t>(phase & 0xFF);
        break;
    }
    uint8_t sample8 = c.wavetable[index];
    if (c.modulation_type == MOD_RING)
    {
        int16_t mod_sample = channels[abs_target_ch].last_sample;
        uint8_t org_sample8 = sample8;
        sample8 = static_cast<uint8_t>CLIP((((static_cast<int>(org_sample8) - 128) * (mod_sample)) * (c.modulation_param_1)) / 65536 + 128, 0, 255);
        //printf("Ring Mod Ch %d Target Ch %d Org %d Sample %d ModSample %d\n", ch, abs_target_ch, org_sample8, sample8, mod_sample);
    }
    else if (c.modulation_type == MOD_SLEW_RATE) {
        int diff = static_cast<int>(sample8 << 8) - static_cast<int>(c.slew_cur);
        c.slew_cur += CLIP(diff, -c.modulation_param_1, c.modulation_param_1); // apply slew rate 
        c.slew_cur = CLIP(c.slew_cur, 0, 255 << 8); // clamp to 0-255
        sample8 = static_cast<uint8_t>(c.slew_cur >> 8);
    }
    else if (c.modulation_type == MOD_SLEW_STEP) {
        int diff = static_cast<int>(sample8 << 8) - static_cast<int>(c.slew_cur);
        c.slew_cur += (diff > 0) ? c.modulation_param_1 : -c.modulation_param_1; // apply slew rate 
        c.slew_cur = CLIP(c.slew_cur, 0, 255 << 8); // clamp to 0-255
        sample8 = static_cast<uint8_t>(c.slew_cur >> 8);
    }
    int16_t out = static_cast<int16_t>((static_cast<int>(sample8) - 128));
    return out;
}

// Every v3 fragment starts with its own predictor state.
void S3W2_Sound::ensureSQV4Frame(Channel& c, uint32_t frame) {
    if (c.sqv4_cached_frame[frame&1]==frame) return;
    const uint8_t* f=pcm_ram->data()+c.pcm_start_addr+16+frame*c.sqv4_header.stride();
    sqv4::State state; state.index=f[1];
    unsigned count=std::min(63u,c.sqv4_header.samples-frame*63);
    for (unsigned i=0;i<count;i++) {
        int sample=sqv4::decodeCode(c.sqv4_header,state,sqv4::codeAt(f+2,i,c.sqv4_header.bits()));
        c.sqv4_frames[frame&1][i]=sqv4::applyVolume(sample,f[0],3);
    }
    c.sqv4_cached_frame[frame&1]=frame;
}

int32_t S3W2_Sound::generateSQV4Sample(int ch)
{
    Channel& c=channels[ch];
    if (!c.sqv4_ready) {
        uint32_t start=c.pcm_start_addr, end=c.pcm_end_addr, loop=c.pcm_loop_addr;
        if (start>=end || end>PCM_RAM_SIZE ||
            !sqv4::inspect(pcm_ram->data()+start,end-start,c.sqv4_header) ||
            c.sqv4_header.channels!=1 || !c.sqv4_header.samples) {
            c.active=false; return 0;
        }
        c.sqv4_loop=0;
        if (loop && loop!=start) {
            if (loop<start+16 || loop>=end || (loop-start-16)%c.sqv4_header.stride()) {
                c.active=false; return 0;
            }
            c.sqv4_loop=((loop-start-16)/c.sqv4_header.stride())*63;
        }
        c.sqv4_state=sqv4::State(); c.sqv4_next_frame=0; c.sqv4_ready=true;
        c.sqv4_cached_frame.fill(UINT32_MAX);
    }
    uint64_t target=static_cast<uint64_t>(c.phase*256);
    if (target>=c.sqv4_header.samples) {
        if (!(c.pcm_control&2)) { c.active=false; return 0; }
        target=c.sqv4_loop+(target-c.sqv4_loop)%(c.sqv4_header.samples-c.sqv4_loop);
        c.phase=(target+(c.phase*256-std::floor(c.phase*256)))/256.0;
    }
    ensureSQV4Frame(c,target/63);
    int32_t current=c.sqv4_frames[(target/63)&1][target%63];
    int32_t q8=current*256;
    int frac=CLIP(static_cast<int>((c.phase*256-std::floor(c.phase*256))*256),0,255);
    if (revisionA && frac && !(c.sample_control&4) && target+1<c.sqv4_header.samples) {
        uint32_t next=target+1;
        ensureSQV4Frame(c,next/63);
        int32_t following=c.sqv4_frames[(next/63)&1][next%63];
        q8=current*(256-frac)+following*frac;
    }
    c.phase+=static_cast<double>(c.frequency)/SOUND_CLOCK/8;
    return q8;
}

int16_t S3W2_Sound::generatePCMSample(int ch)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    Channel &c = channels[ch];
    if (c.waveform_type != WAVEFORM_PCM)
        return 0;
    // bounds check
    if (!pcm_ram)
        return 0;
    uint64_t phase = static_cast<uint64_t>(c.phase*256) & 0xFFFFFFFFFFFFFFFF;
    uint8_t abs_target_ch = convertToAbsoluteChannelAddress(static_cast<uint8_t>(ch), c.modulation_targeting_mode, c.modulation_target);
    if (c.modulation_type == MOD_PHASE)
    {
        int16_t phase_offset = static_cast<int16_t>((c.modulation_param_1 * channels[abs_target_ch].last_sample) >> 12);
        phase = (phase + phase_offset) & 0xFFFFFFFFFFFFFFFF;
    }
    uint32_t addr = (phase+c.pcm_start_addr) & (PCM_RAM_SIZE - 1);
    uint8_t sample8 = (*pcm_ram)[addr];
    if (c.modulation_type == MOD_RING)
    {
        int16_t mod_sample = channels[abs_target_ch].last_sample;
        uint8_t org_sample8 = sample8;
        sample8 = static_cast<uint8_t>CLIP((((static_cast<int>(org_sample8) - 128) * (mod_sample)) * (c.modulation_param_1)) / 65536 + 128, 0, 255);
        //printf("Ring Mod Ch %d Target Ch %d Org %d Sample %d ModSample %d\n", ch, abs_target_ch, org_sample8, sample8, mod_sample);
    }
    else if (c.modulation_type == MOD_SLEW_RATE) {
        int diff = static_cast<int>(sample8 << 8) - static_cast<int>(c.slew_cur);
        c.slew_cur += CLIP(diff, -c.modulation_param_1, c.modulation_param_1); // apply slew rate 
        c.slew_cur = CLIP(c.slew_cur, 0, 255 << 8); // clamp to 0-255
        sample8 = static_cast<uint8_t>(c.slew_cur >> 8);
    }
    else if (c.modulation_type == MOD_SLEW_STEP) {
        int diff = static_cast<int>(sample8 << 8) - static_cast<int>(c.slew_cur);
        c.slew_cur += (diff > 0) ? c.modulation_param_1 : -c.modulation_param_1; // apply slew rate 
        c.slew_cur = CLIP(c.slew_cur, 0, 255 << 8); // clamp to 0-255
        sample8 = static_cast<uint8_t>(c.slew_cur >> 8);
    }
    int16_t out = static_cast<int16_t>((static_cast<int>(sample8) - 128));
    // advance
    if (addr >= c.pcm_end_addr)
    {
        if (c.pcm_control & 0x02)
        { // loop on 
            c.phase = (c.pcm_loop_addr-c.pcm_start_addr)/256.0;
        }
        else
        {
            c.phase = (c.pcm_end_addr-c.pcm_start_addr)/256.0;
            out = 0; // silence after end
        }
    } else {
        c.phase+=static_cast<double>(c.frequency)/SOUND_CLOCK/8; // advance phase
    }
    //printf("PCM Ch %d curAddr %05X startAddr %05X endAddr %05X loopAddr %05X Sample %d\n", ch, addr, c.pcm_start_addr, c.pcm_end_addr, c.pcm_loop_addr, out);
    return out;
}

#define MOD_NOISE_CUSTOM_TAP MOD_PHASE

int16_t S3W2_Sound::generateNoiseSample(int ch)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    Channel &c = channels[ch];
    if (c.waveform_type != WAVEFORM_NOISE)
        return 0;
    // simple 1-bit LFSR -> -32768 or +32767
    // advance
    c.phase+=static_cast<double>(c.frequency)/SOUND_CLOCK; // advance phase
    uint64_t phase = static_cast<uint64_t>(c.phase*256) & 0xFFFFFFFFFFFFFFFF;
    if (c.old_phase/8 != phase/8) {
        c.old_phase = static_cast<uint64_t>(c.phase*256) & 0xFFFFFFFFFFFFFFFF;
        uint32_t bit = 0;
        uint32_t taps = 0;
        switch (c.modulation_type)
        {
            case MOD_NONE:
                bit = ((c.lfsr_state >> 0) ^ (c.lfsr_state >> 1)) & 1;
                c.lfsr_state = (c.lfsr_state >> 1) | (bit << 22); // update LFSR
                break;
            case MOD_NOISE_CUSTOM_TAP:
                taps = c.modulation_param_1 << 16 | c.modulation_param_2; // as 23-bit tap mask (ignore upper bits)
                bit = 0;
                for (int i = 0; i < 23; ++i)
                {
                    if (taps & (1 << i))
                        bit ^= (c.lfsr_state >> i) & 1;
                }
                c.lfsr_state = (c.lfsr_state >> 1) | (bit << 22); // update LFSR
                break;
            default:
                bit = ((c.lfsr_state >> 0) ^ (c.lfsr_state >> 1)) & 1;
                c.lfsr_state = (c.lfsr_state >> 1) | (bit << 22); // update LFSR
                break;
        }
    }
    uint8_t abs_target_ch = convertToAbsoluteChannelAddress(static_cast<uint8_t>(ch), c.modulation_targeting_mode, c.modulation_target);
    uint8_t sample8 = (c.lfsr_state & 1) ? 255 : 0;
    if (c.modulation_type == MOD_RING)
    {
        int16_t mod_sample = channels[abs_target_ch].last_sample;
        uint8_t org_sample8 = sample8;
        sample8 = static_cast<uint8_t>CLIP((((static_cast<int>(org_sample8) - 128) * (mod_sample)) * (c.modulation_param_1)) / 65536 + 128, 0, 255);
        //printf("Ring Mod Ch %d Target Ch %d Org %d Sample %d ModSample %d\n", ch, abs_target_ch, org_sample8, sample8, mod_sample);
    }
    else if (c.modulation_type == MOD_SLEW_RATE) {
        int diff = static_cast<int>(sample8 << 8) - static_cast<int>(c.slew_cur);
        c.slew_cur += CLIP(diff, -c.modulation_param_1, c.modulation_param_1); // apply slew rate 
        c.slew_cur = CLIP(c.slew_cur, 0, 255 << 8); // clamp to 0-255
        sample8 = static_cast<uint8_t>(c.slew_cur >> 8);
    }
    else if (c.modulation_type == MOD_SLEW_STEP) {
        int diff = static_cast<int>(sample8 << 8) - static_cast<int>(c.slew_cur);
        c.slew_cur += (diff > 0) ? c.modulation_param_1 : -c.modulation_param_1; // apply slew rate 
        c.slew_cur = CLIP(c.slew_cur, 0, 255 << 8); // clamp to 0-255
        sample8 = static_cast<uint8_t>(c.slew_cur >> 8);
    }
    int16_t out = static_cast<int16_t>((static_cast<int>(sample8) - 128));
    return out;
}

int16_t S3W2_Sound::generateDMAPCMSample(int ch)
{
    // not implemented yet, raise error
    throw std::runtime_error("DMA PCM not implemented");
}

int16_t S3W2_Sound::applyModulation(int ch, int16_t &sample)
{
    // placeholder: no modulation applied yet
    return static_cast<int16_t>(sample);
}

void S3W2_Sound::mixChannels(int16_t *left, int16_t *right)
{
    // placeholder
    (void)left;
    (void)right;
}

void S3W2_Sound::writeChannelWavetable(int ch, uint8_t offset, uint8_t value)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return;
    if (offset >= WAVETABLE_SIZE)
        return;
    channels[ch].wavetable[offset] = value;
}

void S3W2_Sound::writeChannelPCMReg(int ch, uint8_t offset, uint8_t value)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return;
    Channel &c = channels[ch];
    if (offset<=8 || (offset==9 && !(c.pcm_control&1) && (value&1))) c.sqv4_ready=false;
    // offset mapping: 0-2 start addr, 3-5 end, 6-8 loop, 9 control
    if (offset <= 2)
    {
        // start addr is 20-bit across 3 bytes big-endian. We'll write per-byte
        uint32_t shift = (2 - offset) * 8;
        uint32_t mask = 0xFFu << shift;
        uint32_t v = (static_cast<uint32_t>(value) << shift);
        c.pcm_start_addr = (c.pcm_start_addr & ~mask) | v;
    }
    else if (offset >= 3 && offset <= 5)
    {
        uint32_t shift = (5 - offset) * 8;
        uint32_t mask = 0xFFu << shift;
        uint32_t v = (static_cast<uint32_t>(value) << shift);
        c.pcm_end_addr = (c.pcm_end_addr & ~mask) | v;
    }
    else if (offset >= 6 && offset <= 8)
    {
        uint32_t shift = (8 - offset) * 8;
        uint32_t mask = 0xFFu << shift;
        uint32_t v = (static_cast<uint32_t>(value) << shift);
        c.pcm_loop_addr = (c.pcm_loop_addr & ~mask) | v;
    }
    else if (offset == 0x0A)
    {
        // 3WS8PNB SAMPLE_CONTROL; bits 7:3 reserved. Phase/SRAM are untouched.
        if (revisionB) c.sample_control = value & 0x07;
    }
    else if (offset == 9)
    {
        c.pcm_control = value;
        if (value & 0x01)
        { // play
            c.active = true;
        }
        else
        {
            c.active = false;
        }
    }
}

void S3W2_Sound::writePCMRAM(uint32_t address, uint8_t value)
{
    if (!pcm_ram)
        return;
    if (address >= PCM_RAM_SIZE)
        return;
    (*pcm_ram)[address] = value;
    for (auto& c: channels) c.sqv4_ready=false;
}

uint8_t S3W2_Sound::readPCMRAM(uint32_t address)
{
    if (!pcm_ram)
        return 0;
    if (address >= PCM_RAM_SIZE)
        return 0;
    return (*pcm_ram)[address];
}

void S3W2_Sound::writeChannelControl(int ch, uint8_t offset, uint8_t value)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return;
    Channel &c = channels[ch];
    if (offset == 0)
    {
        // frequency MSB (big endian across 0-1)
        c.frequency = (c.frequency & 0x00FF) | (static_cast<uint16_t>(value) << 8);
    }
    else if (offset == 1)
    {
        c.frequency = (c.frequency & 0xFF00) | static_cast<uint16_t>(value);
    }
    else if (offset == 2)
    {
        if (!revisionA && value==WAVEFORM_FIFO_DAC) value=WAVEFORM_WAVETABLE;
        if (c.waveform_type!=value) c.sqv4_ready=false;
        c.waveform_type = value;
    }
    else if (offset == 3)
    {
        c.volume = value;
    }
    else if (offset == 4)
    {
        c.panpot = value;
    }
    else if (offset == 5)
    {
        c.modulation_type = static_cast<uint8_t>(value >> 3);
        c.modulation_target = (value & 0x07);
    }
    else if (offset >= 6 && offset <= 9)
    {
        if (offset<=7) c.modulation_param_1 = (c.modulation_param_1 & ~(0xFFu << ((7 - offset) * 8))) |
                             (static_cast<uint16_t>(value) << ((7 - offset) * 8));
        else c.modulation_param_2 = (c.modulation_param_2 & ~(0xFFu << ((9 - offset) * 8))) |
                             (static_cast<uint16_t>(value) << ((9 - offset) * 8));
    }
    else if (offset == 0x0A)
    {
        // access resets phase
        c.phase = 0.0;
        c.sqv4_ready = false;
        c.slew_cur = 0;
        c.lfsr_state = 0x12D4803C;
        resetFIFO(c);
    }
    else if (offset == 0x0B)
    {
        c.modulation_targeting_mode = value & 0x01;
    }
    else if (revisionA && offset==0x0C) c.fifo_latch=value;
    else if (revisionA && offset==0x0D && c.fifo.size()<2048) {
        c.fifo.push_back(static_cast<int16_t>((static_cast<uint16_t>(c.fifo_latch)<<8)|value));
    }
}

uint8_t S3W2_Sound::readChannelWavetable(int ch, uint8_t offset)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    if (offset >= WAVETABLE_SIZE)
        return 0;
    return channels[ch].wavetable[offset];
}

uint8_t S3W2_Sound::readChannelPCMReg(int ch, uint8_t offset)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    Channel &c = channels[ch];
    if (offset <= 2)
    {
        uint32_t shift = (2 - offset) * 8;
        return static_cast<uint8_t>((c.pcm_start_addr >> shift) & 0xFF);
    }
    else if (offset >= 3 && offset <= 5)
    {
        uint32_t shift = (5 - offset) * 8;
        return static_cast<uint8_t>((c.pcm_end_addr >> shift) & 0xFF);
    }
    else if (offset >= 6 && offset <= 8)
    {
        uint32_t shift = (8 - offset) * 8;
        return static_cast<uint8_t>((c.pcm_loop_addr >> shift) & 0xFF);
    }
    else if (offset == 9)
    {
        return c.pcm_control;
    }
    return 0;
}

uint8_t S3W2_Sound::readChannelControl(int ch, uint8_t offset)
{
    if (ch < 0 || ch >= NUM_CHANNELS)
        return 0;
    Channel &c = channels[ch];
    if (offset == 0)
    {
        return static_cast<uint8_t>((c.frequency >> 8) & 0xFF);
    }
    else if (offset == 1)
    {
        return static_cast<uint8_t>(c.frequency & 0xFF);
    }
    else if (offset == 2)
    {
        return c.waveform_type;
    }
    else if (offset == 3)
    {
        return c.volume;
    }
    else if (offset == 4)
    {
        return c.panpot;
    }
    else if (offset == 5)
    {
        return static_cast<uint8_t>((c.modulation_type << 3) | (c.modulation_target & 0x07));
    }
    else if (offset >= 6 && offset <= 9)
    {
        uint32_t shift = ((offset<=7?7:9) - offset) * 8;
        if (offset <= 7)
        {
            return static_cast<uint8_t>((c.modulation_param_1 >> shift) & 0xFF);
        }
        else
        {
            return static_cast<uint8_t>((c.modulation_param_2 >> shift) & 0xFF);
        }
    }
    else if (offset == 0x0A)
    {
        // access resets phase
        c.phase = 0.0;
        c.sqv4_ready = false;
        return 0;
    }
    else if (offset == 0x0B)
    {
        return c.modulation_targeting_mode & 0x01;
    }
    else if (revisionA && (offset==0x0E || offset==0x0F)) {
        return static_cast<uint8_t>(c.fifo.size()>>((0x0F-offset)*8));
    }
    else if (offset>=0x10 && offset<=0x19) return readChannelPCMReg(ch,offset-0x10);
    else if (revisionB && offset==0x1A) return c.sample_control;
    return 0;
}

uint32_t S3W2_Sound::extractAddress20bit(uint32_t reg_addr)
{
    // reg_addr assumed to point to first of three bytes; big endian
    // Not implemented in this header-only context; kept for API compatibility
    (void)reg_addr;
    return 0;
}

void S3W2_Sound::setAddress20bit(uint32_t &target, uint32_t reg_addr, uint8_t value)
{
    (void)reg_addr;
    (void)value;
    (void)target;
}

// Convert modulation target to absolute channel address
uint8_t S3W2_Sound::convertToAbsoluteChannelAddress(uint8_t carrier_channel, uint8_t modulation_targeting_mode, uint8_t modulation_target)
{
    if (modulation_targeting_mode == 0)
    {
        // absolute addressing
        return modulation_target & 0x07;
    }
    else if (modulation_targeting_mode == 1)
    {
        // relative addressing
        int8_t relative_offset = static_cast<int8_t>(modulation_target & 0x07);
        if (relative_offset >= 4)
            relative_offset -= 8; // convert to signed 3-bit
        int absolute_channel = static_cast<int>(carrier_channel) + relative_offset;
        if (absolute_channel < 0)
            absolute_channel += NUM_CHANNELS;
        else if (absolute_channel >= NUM_CHANNELS)
            absolute_channel -= NUM_CHANNELS;
        absolute_channel &= 0x07; // ensure within 0-7 (wraps around)
        return static_cast<uint8_t>(absolute_channel);
    }
    return modulation_target&7;
}
