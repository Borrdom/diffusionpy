import matplotlib
import matplotlib.pyplot as plt

def move_figure(f, x, y):
    """Move figure's upper left corner to pixel (x, y)"""
    # check if 
    backend = matplotlib.get_backend()
    if backend == 'TkAgg':
        f.canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
    elif backend == 'WXAgg':
        f.canvas.manager.window.SetPosition((x, y))
    elif backend=="agg":
        None
    else:
        # This works for QT and GTK
        # You can also use window.setGeometry
        f.canvas.manager.window.move(x, y)

if __name__=="__main__":
    f, ax = plt.subplots()
    move_figure(f, 500, 500)    
    plt.show()