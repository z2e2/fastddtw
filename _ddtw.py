import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def diff_dist(x,y):
    '''
    Computing drivative differences between dx and dy
    Arguments:
        x -- dx from signal 1, numpy array of shape ( 3,  )
        y -- dy from signal 2, numpy array of shape ( 3,  )
    Result:
          -- absolute difference of estimated derivatives of x, y
    '''
    dx = ((x[1] - x[0]) + (x[2]-x[0])/2)/2
    dy = ((y[1] - y[0]) + (y[2]-y[0])/2)/2
    return abs(dx-dy)

def DDTW(signal_1, signal_2):
    '''
    Arguments:
        signal_1 -- first time series, numpy array of shape ( n1,  )
        signal_2 -- second time series, numpy array of shape ( n2,  )
    Results:
        ddtw -- distance matrix, numpy array of shape ( n1 - 2, n2 - 2 )
        ddtw_traceback -- traceback matrix, numpy array of shape ( n1 - 2, n2 - 2 )
    ''' 
    assert signal_1.shape[0] != 0 and signal_2.shape[0] != 0, '''Input signals must be a column vectors,
                                                                Please check the input signal dimension.'''
    assert signal_1.shape[0] >= 3 and signal_2.shape[0] >= 3, '''The length of your signal should be 
                                                                 greater than 3 to implement DDTW.'''
    n_rows = signal_1.shape[0]-2
    n_cols = signal_2.shape[0]-2
    ddtw = np.zeros((n_rows,n_cols))
    ddtw_traceback = np.zeros((n_rows,n_cols))
    ddtw[0,0] = diff_dist(signal_1[0:3], signal_2[0:3])
    for i in range(1, n_rows):
        ddtw[i,0] = ddtw[i-1,0] + diff_dist(signal_1[i-1:i+2], signal_2[0:3])
        ddtw_traceback[i,0] = 1
    for j in range(1, n_cols):
        ddtw[0,j] = ddtw[0,j-1] + diff_dist(signal_1[0:3], signal_2[j-1:j+2])
        ddtw_traceback[0,j] = 2
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            temp = np.array([ddtw[i-1,j-1], ddtw[i-1,j], ddtw[i,j-1]])
            best_idx = np.argmin(temp)
            ddtw[i,j] = diff_dist(signal_1[i-1:i+2], signal_2[j-1:j+2]) + temp[best_idx]
            ddtw_traceback[i,j] = best_idx
    return ddtw, ddtw_traceback

def plot_traceback(ddtw,ddtw_traceback):
    i, j = ddtw.shape
    i -= 1
    j -= 1
    x = [i]
    y = [j]
    while (i != 0 or j != 0):
        if i != 0 and j != 0:
            idx_i = [i-1, i-1, i]
            idx_j = [j-1 ,j, j-1]
            idx = int(ddtw_traceback[i,j])
            i = idx_i[idx]
            j = idx_j[idx]
        elif i == 0 and j != 0:
            j = j-1
        elif i != 0 and j == 0:
            i = i-1
        elif i==1 and j == 1:
            i = 0
            j = 0
        x.append(i)
        y.append(j)
    plt.plot(x, y)
    plt.show()
    return x, y
def plot_aligned_sig(y1, y2, x, y):
    new_x = np.zeros(len(x))
    new_y = np.zeros(len(y))
    for i in range(len(x)):
        new_x[i] = y1[x[i]+1]
        new_y[i] = y2[y[i]+1]
    plt.plot(new_x[::-1])
    plt.plot(new_y[::-1])
    plt.title('signal 1')
    plt.show()
    return new_x[::-1], new_y[::-1]
def plot_raw_signals(signal_1, signal_2, title = 'raw_signals'):
    '''
    Arguments:
        signal_1 -- first time series, numpy array of shape ( n1,  )
        signal_2 -- second time series, numpy array of shape ( n2,  )
    Results:
          Figure 
    ''' 
    plt.plot(signal_1)
    plt.plot(signal_2)
    plt.grid()
    plt.title(title)
    plt.xlabel('time')
    plt.ylabel('value')
    plt.show()