# register view

during playback, "Register View" shows the chips' registers (data involved with each operation).

values that change between display frames flash light red (#FF4F4F) and fade back to the normal text color over half a second. another change restarts the highlight. writes that leave the displayed value unchanged do not trigger a highlight.

![register view window](register.png)

the menu contains the following options:
- **Bytes per column**: sets the number of columns displayed. useful for aligning registers for convenient viewing and comparison.
- **Highlight register changes**: enables or disables the fading highlight for changed values.
- **Highlight color**: changes the highlight color, including by hexadecimal RGB value. defaults to #FF4F4F. these settings are saved between sessions.
