#!/usr/bin/env python
"""
 Blinks an LED on digital pin 13
 in 1 second intervals
"""

from Arduino import Arduino
import time
import matplotlib.pyplot as plt

board = Arduino(baud=9600,port="COM6") # plugged in via USB, serial com at rate 115200
board.pinMode(13, "OUTPUT")
board.pinMode(2, "INPUT")

i=0
while True:
    i+=1
    board.digitalWrite(13, "LOW")
    time.sleep(1)
    board.digitalWrite(13, "HIGH")
    time.sleep(1)
    plt.plot(i,board.analogRead(14),'kx')
    plt.pause(0.05)
plt.show()