import numpy as np

def dendShape( saddle_start_um, theta, rm, rp ):
         # total dome length = width = 2 * Rd * sin(theta), for theta = 0.01
         dth = theta / 100
         xbend = 2 * np.sin( theta ) * (rp + rm )

         Lflat = saddle_start_um * 1e6

         x = [0.0, Lflat]
         y = [0.0, 0.0]
         L_truncate = 2 * Lflat + xbend
    
    #Starting from flat through the saddle
         x.extend( [ Lflat + rp * np.sin( th + dth ) for th in np.arange(0.0, theta * 0.999999, dth ) ] )
         y.extend( [ rp * (1-np.cos( th + dth) ) for th in np.arange(0.0, theta * 0.999999, dth ) ] )
         dome_start = x[-1] * 1e-6
    #End of the saddle
         xlast = x[-1]
         ylast = y[-1]
    #Going to the middle of dome
         xoffset = rm * np.sin( theta ) + xlast
         yoffset = -rm * np.cos( theta ) + ylast

    #Going from the beginning of the dome to the middle.
         x.extend( [ xoffset - rm * np.sin( th ) for th in np.arange (theta, 0, -dth ) ] )
         y.extend ([yoffset + rm * np.cos( th ) for th in np.arange (theta, 0, -dth ) ] )
         xlast = x[-1]
         ylast = y[-1]
         x.extend( [ L_truncate - i for i in x[::-1] ] )
         y.extend( y[::-1] )
         dome_end = ( xoffset + rm * np.sin( theta ) ) * 1e-6
         #dome_end = x[-1]
         xlast = x[-1]
         ylast = y[-1]

         return np.array( x ), np.array( y ), dome_start, dome_end

def get_shape(saddle_start_list, theta_list, Hmean_list, rp):
         x2 = np.asarray([])
         y2 = np.asarray([])
         for sh in range(len(saddle_start_list)):
             x, y, dome_start, dome_end = dendShape( saddle_start_list[sh] , theta_list[sh], ( 1 / Hmean_list[sh] ) * 1e6, rp[sh] * 1e6 )
             x2 = np.append(x2, x)
             y2 = np.append(y2, y)
             x2_trunc = np.where(x2 <= 20)
         return x2[x2_trunc], y2[x2_trunc]



