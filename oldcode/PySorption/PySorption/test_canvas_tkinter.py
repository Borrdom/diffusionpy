# import tkinter
# import matplotlib.pyplot as plt
# import numpy as np
# import matplotlib
# from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# matplotlib.use("TkAgg")
# x=np.linspace(0,10,50)
# y=x**2
# fig,ax=plt.subplots()
# ax.plot(x,y)
# top = tkinter.Tk()
# canvas = FigureCanvasTkAgg(fig, master=top)
# canvas.get_tk_widget().pack()
# canvas.draw()
# top.mainloop()

import tkinter
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

class App:
    def __init__(self, master):
        # Create a container
        frame = tkinter.Frame(master)
        # Create 2 buttons
        self.button_left = tkinter.Button(frame,text="< Decrease Slope",
                                        command=self.decrease)
        self.button_left.pack(side="left")
        self.button_right = tkinter.Button(frame,text="Increase Slope >",
                                        command=self.increase)
        self.button_right.pack(side="left")

        fig = Figure()
        ax = fig.add_subplot(111)
        self.line, = ax.plot(range(10))

        self.canvas = FigureCanvasTkAgg(fig,master=master)
        #self.canvas.show()
        self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)
        frame.pack()

    def decrease(self):
        x, y = self.line.get_data()
        self.line.set_ydata(y - 0.2 * x)
        self.canvas.draw()

    def increase(self):
        x, y = self.line.get_data()
        self.line.set_ydata(y + 0.2 * x)
        self.canvas.draw()

root = tkinter.Tk()
app = App(root)
root.mainloop()