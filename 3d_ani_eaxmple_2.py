
def animateGraph(self):
        global fig
        fig = plt.figure()
        ax1 = fig.add_subplot(1,1,1)
        # animate()
        # Description: Every 200 ms, get speed, steering angle, and displacement estimate and update dynamic graph
        def animate(i):
            lx,ly = self.getLocation()
            try:
                ax1.clear()
                ax1.plot(lx,ly)
                ax1.set_title("2D position estimate")
                ax1.set_ylabel(" Y displacement (m)")
                ax1.set_xlabel(" X displacement (m)")
            except:
                print('s')
        plt.grid(True)
        plt.subplots_adjust(hspace = 1,wspace = 0.6)
        ani = animation.FuncAnimation(fig, animate, interval=200)
        plt.show() 