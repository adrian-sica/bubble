import numpy as np
import matplotlib.pyplot as plt
from fluidfoam import readmesh
from fluidfoam import readvector, readscalar

def main():
    
    plt.style.use('mystyle.mplstyle')
    
    n_x = [80, 80, 160]
    max_co = [0.1, 0.01, 0.01]
    case_names = []
    colors = ['black', 'blue', 'green', 'red']
    lineStyles = ['-', '--', ':', '-.']
    
    lx = 1
    ly = 2
    lz = 0.1    
    d_bubble = 0.5
    circle_bubble_area = np.pi*d_bubble*lz
    
    t_start = 0
    t_end = 3
    wrtite_time = 0.05
    times = t_start + np.arange(0, t_end/wrtite_time+1)*wrtite_time
    
    y_label = ['Bubble area (m$^2$)','Bubble circularity',
               'Bubble center of mass (m)','Bubble rise velocity (m/s)']
    fig_name = ['area','circularity','centre','velocity','shape']
    
    case_dir = []
    benchmark = 'c2g2l3'
    benchmark_data = np.loadtxt('../Benchmark/data_bench_quantities/' + benchmark + '.txt')
    benchmark_shape = np.loadtxt('../Benchmark/data_bubble_shapes/' + benchmark + 's.txt')
    
    bubble_area = np.pi*d_bubble**2/4*np.ones((times.size, len(n_x)))
    bubble_circularity = np.ones((times.size, len(n_x)))
    bubble_centre = 0.5*np.ones((times.size, len(n_x)))
    bubble_velocity = np.zeros((times.size, len(n_x)))
        
    for i in range(len(n_x)):
        cell_area = lx*ly/(2*n_x[i]**2)
        case_dir.append('../mesh_' + str(n_x[i]) + '_Co_' + str(max_co[i]))
        _, y, _ = readmesh(case_dir[i])
        
        for j in range(len(times)):
            time_name = (str(np.round(times[j],8)) + '/').replace('.0/','/')
            u = readvector(case_dir[i], time_name, 'U')
            alpha_water = readscalar(case_dir[i], time_name, 'alpha.water')
            alpha_air = 1 - alpha_water
            
            sum_alpha = np.sum(alpha_air)
            sum_y_alpha = np.sum(np.multiply(y, alpha_air))
            sum_uy_alpha = np.sum(np.multiply(u[1], alpha_air))
            
            bubble_area[j, i] = sum_alpha*cell_area
            bubble_centre[j, i] = sum_y_alpha / sum_alpha
            bubble_velocity[j, i] = sum_uy_alpha / sum_alpha
            
        area = np.loadtxt(case_dir[i] + '/postProcessing/bubbleArea/0/'
                          +'surfaceFieldValue.dat', skiprows=5)
        bubble_circularity[1:, i] = np.divide(circle_bubble_area, area[:,1])

    results = [bubble_area, bubble_circularity, bubble_centre, bubble_velocity]
    
    # Plots
    for i in range(5):
        if i < 4: # bubble charactristics
            fig, ax = plt.subplots()
            for j in range(len(n_x)):
                case_names.append('Mesh $' + str(n_x[j]) + r'\times ' + str(2*n_x[j]) 
                                  + r'$, \verb|maxCo|=$' + str(max_co[j]) + '$')
                ax.plot(times, results[i][:,j], linestyle=lineStyles[j], color=colors[j],
                            label=case_names[j])
            ax.plot(benchmark_data[1::15,0], benchmark_data[1::15,i+1], marker='o', 
                    linestyle='None', color='r', markerfacecolor='None', label='FreeLIFE')
            ax.set_xlabel('Time (s)')
            ax.set_ylabel(y_label[i])
            ax.set_xlim([0, 3])
        else: # bubble shapes
            fig, ax = plt.subplots(figsize=(6.2,7))
            for j in range(len(n_x)):
                shape = np.loadtxt(case_dir[j] + '/postProcessing/interface/3/'
                                    +'alpha.water_interface.raw', skiprows=2)
                ax.plot(shape[:,0], shape[:,1], linestyle='None', color=colors[j], marker='.',
                        markersize=3, label=case_names[j])
                
            ax.plot(benchmark_shape[1::15,0], benchmark_shape[1::15,1], linestyle='None',
                    color='r', marker='.', markersize=3, label='FreeLIFE')
            ax.set_xlabel('$x$ (m)')
            ax.set_ylabel('$y$ (m)')
            ax.axis('equal')
        ax.legend(loc='best', edgecolor='k', labelspacing = 0.25)
        fig.savefig('Figures/' + fig_name[i] + '.pdf', bbox_inches='tight')
        
    plt.show()

if __name__ == '__main__':
    main()
