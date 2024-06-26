#include <math.h>
#include <random>


class Cpt100_sound {
public:
    
    #include "envelove.cpp"
    #include "ram.cpp"
    std::random_device rd;
    std::mt19937 mt;
    long long Total_time = 0;
    double t1[4] = {0,0,0,0};
    double t2[4] = {0,0,0,0};
    double t3[4] = {0,0,0,0};
    double t4[4] = {0,0,0,0};
    double twt[2] = {0,0};
    double in1[2]  = {0.0,0.0};
    double in2[2]  = {0.0,0.0};
    double out1[2] = {0.0,0.0};
    double out2[2] = {0.0,0.0};
    std::vector<int> gateTick = {0,0,0,0};
    std::vector<Byte> reg,regenvl,regwt;
    std::vector<int> noise;
    std::vector<std::vector<double>> sintable;
    std::vector<EnvGenerator> envl;
    EnvGenerator _envl;
    double prev = 0;
    int pcm_addr[2],pcm_len[2],pcm_loop[2]={0,0};
    std::vector<std::vector<Byte>> pcm_ram; 

    Cpt100_sound() {
    };
    
    #define CPT100_SAMPLE_FREQ 48000

    double sind(double theta) {
        return sin(theta*2*M_PI);
    }

    double modulate(double theta, int wf) {
        return sintable.at(wf).at((int)((theta*256))&0xff);
    }

    double generateFMWave(double t1, double v1, double t2, double v2, double t3, double v3, double t4, double v4, int w1, int w2, int w3, int w4) {

        double value = modulate(t1+modulate(t2+modulate(t3+modulate(t4,w4)*v4,w3)*v3,w2)*v2,w1)*v1*255*127;
        return value;

    }

    void applyEnveloveToRegisters(std::vector<Byte> &reg, std::vector<Byte> &regenvl, int opNum, int ch, double dt) {
        ADSRConfig adsr;
        adsr.attackTime = ((double)regenvl.at(ch*16+opNum*4+0).toInt())/64;
        adsr.decayTime = ((double)regenvl.at(ch*16+opNum*4+1).toInt())/64;
        adsr.sustainLevel = ((double)regenvl.at(ch*16+opNum*4+2).toInt())/255;
        adsr.releaseTime = ((double)regenvl.at(ch*16+opNum*4+3).toInt())/64;
        if (regenvl.at(64+ch).toInt() == 0 && gateTick.at(ch) == 1) {
            envl.at((size_t)(ch*4+opNum)).noteOff();
            if(opNum == 3) {
                gateTick.at(ch)=0;
            }
            
        }
        if (regenvl.at(64+ch).toInt() == 1 && gateTick.at(ch) == 0) {
            envl.at((size_t)(ch*4+opNum)).reset(EnvGenerator::State::Attack); 
            if(opNum == 3) {
                gateTick.at(ch)=1;
            }
        }
        //std::cout << dt << std::endl; //envl.at((size_t)(ch*4+opNum)).m_elapsed
        reg.at(ch*16+opNum+5) = (Byte)(envl.at((size_t)(ch*4+opNum)).currentLevel()*255*((double)(reg.at(ch*16+opNum+9).toInt())/255));
        envl.at((size_t)(ch*4+opNum)).update(adsr,dt);
    }

    std::vector<std::vector<int16_t>> AudioCallBack(int len)
    {
        int i;
        std::vector<int16_t> _frames(len,0);
        std::vector<std::vector<int16_t>> frames(6,_frames);
        int framesize = len;
        reg = ram_peek2array(ram,0x10000,64);
        regenvl = ram_peek2array(ram,0x10040,68);
        regwt = ram_peek2array(ram,0x10084,76);

        for (int ch=0;ch<2;ch++) {
            if(regwt.at(ch+6).toInt() == 4) {
                pcm_addr[ch] = regwt.at(12+32*ch+0).toInt()*65536+regwt.at(12+32*ch+1).toInt()*256+regwt.at(12+32*ch+2).toInt();
                pcm_len[ch] = regwt.at(12+32*ch+3).toInt()*65536+regwt.at(12+32*ch+4).toInt()*256+regwt.at(12+32*ch+5).toInt();
                pcm_loop[ch] = regwt.at(12+32*ch+6).toInt()*65536+regwt.at(12+32*ch+7).toInt()*256+regwt.at(12+32*ch+8).toInt();
                //std::cout << pcm_addr[ch] << std::endl;
                //std::cout << pcm_len[ch] << std::endl;
            }
        }
        
        for(int ch=0; ch < 4; ch++) {
            for (int opNum=0; opNum < 4; opNum++) {
                applyEnveloveToRegisters(reg,regenvl,opNum,ch,((double)len/(double)CPT100_SAMPLE_FREQ));
                ram_poke(ram,0x10000+16*ch+opNum+5,reg.at(16*ch+opNum+5));
            }
        }
        for (i = 0; i < framesize; i++) {
            double result[6] = {0};
            for(int ch=0; ch < 4; ch++) {
                int addr = 16*ch;
                double f1 = ((double)reg.at(addr+0).toInt()*256+reg.at(addr+1).toInt());
                t1[ch] = t1[ch] + (f1/CPT100_SAMPLE_FREQ);
                double v1 = ((double)reg.at(addr+5).toInt())/255;
                t2[ch] = t2[ch] + ((double)(f1*reg.at(addr+2).toInt()))/16/CPT100_SAMPLE_FREQ;
                double v2 = ((double)reg.at(addr+6).toInt())/128;
                t3[ch] = t3[ch] + ((double)(f1*reg.at(addr+3).toInt()))/16/CPT100_SAMPLE_FREQ;
                double v3 = ((double)reg.at(addr+7).toInt())/128;
                t4[ch] = t4[ch] + ((double)(f1*reg.at(addr+4).toInt()))/16/CPT100_SAMPLE_FREQ;
                double v4 = ((double)reg.at(addr+8).toInt())/128;
                int w1 = reg.at(addr+13).toInt()>>4;
                int w2 = reg.at(addr+13).toInt()&0xf;
                int w3 = reg.at(addr+14).toInt()>>4;
                int w4 = reg.at(addr+14).toInt()&0xf;
                result[ch] += generateFMWave(t1[ch],v1,t2[ch],v2,t3[ch],v3,t4[ch],v4,w1,w2,w3,w4);
            }
            for(int ch=0; ch<2; ch++) {
                double ft = ((double)regwt.at(ch*2+0).toInt()*256+regwt.at(ch*2+1).toInt());
                twt[ch] = twt[ch] + (ft/CPT100_SAMPLE_FREQ)*32;
                double vt = ((double)regwt.at(ch+4).toInt())/255;
                int val = 0;
                if (regwt.at(ch+6).toInt() == 1) {
                    val = regwt.at(12+32*ch+((int)twt[ch]%32)).toInt();
                    if ((int)(twt[ch]*2)%2 == 0) {
                        val /= 16;
                    } else {
                        val %= 16;
                    }
                    val *= 16;
                } else if(regwt.at(ch+6).toInt() == 0) {
                val = regwt.at(12+32*ch+((int)twt[ch]%32)).toInt();
                } else if(regwt.at(ch+6).toInt() == 2) {
                    val = noise.at(((int)twt[ch]%65536))*255;
                } else if(regwt.at(ch+6).toInt() == 3) {
                    val = noise.at(((int)twt[ch]%64))*255;
                } else if(regwt.at(ch+6).toInt() == 4) {
                    if (pcm_addr[ch]+(int)twt[ch] > pcm_len[ch] && pcm_loop[ch] < pcm_len[ch]) {
                        val = ram_peek(ram,pcm_addr[ch]+((int)twt[ch]%(pcm_len[ch]-pcm_loop[ch]))+(pcm_addr[ch]-pcm_loop[ch])).toInt();
                    } else {
                        val = ram_peek(ram,std::min(pcm_addr[ch]+(int)twt[ch],pcm_len[ch])).toInt();
                    }
                }
                val -= 128;
                double omega = 2.0 * 3.14159265 * ((double)regwt.at(ch+8).toInt()+1)*32 / CPT100_SAMPLE_FREQ;
                double alpha = sin(omega) / (2.0 * 1.5);
                double a0 =  1.0 + alpha;
                double a1 = -2.0 * cos(omega);
                double a2 =  1.0 - alpha;
                double b0 = (1.0 - cos(omega)) / 2.0;
                double b1 =  1.0 - cos(omega);
                double b2 = (1.0 - cos(omega)) / 2.0;
                double output = b0/a0*(double)val+b1/a0*in1[ch]+b2/a0*in2[ch]-a1/a0*out1[ch]-a2/a0*out2[ch];
                in2[ch]  = in1[ch];
                in1[ch]  = val; 
                out2[ch] = out1[ch];     
                out1[ch] = output; 
                if (regwt.at(ch+8).toInt() == 0) {
                    result[ch+4] += (double)(val)*255*vt;
                } else {
                    result[ch+4] += std::min(std::max((double)output*255*vt,-32768.0),32767.0);
                }
            }
            for(int ch=0; ch<6; ch++) {
                frames[ch][i] = (int16_t)result[ch];
            }
        }

        reg.clear();
        regenvl.clear();
        regwt.clear();
        Total_time++;
        return frames;
    }

    void initSound() {
        
        envl.resize(16,_envl);
        noise.resize(65536,0);
        std::vector<double> _sintable;
        _sintable.resize(256,0);
        sintable.resize(16,_sintable);
        for (int i=0;i<65536;i++) {
            noise[i] = mt()%2;
        }
        for (int i=0; i<256; i++) {
            sintable.at(0).at(i) = sind((double)i/256);
            sintable.at(1).at(i) = std::max(0.0,sind((double)i/256));
            sintable.at(4).at(i) = i<128?(double)i/64-1.0:(double)i/-64+3.0;
            sintable.at(5).at(i) = (double)i/128-1.0;
            if (i<128) {
                sintable.at(2).at(i) = sind((double)i/128);
                sintable.at(3).at(i) = 1.0;
            } else {
                sintable.at(3).at(i) = -1.0;
            }
        }
        
        for (int addr=0x10090;addr<0x100d0;addr++) {
            ram_poke(ram,addr,0x00);
        }
    }

    void resetGate(int ch) {
        for (int i=0;i<4;i++) {
            envl.at((size_t)(ch*4+i)).reset(EnvGenerator::State::Attack); 
        }
    }

    void wtSync(int ch) {
        twt[ch]=0;
    }
    
};