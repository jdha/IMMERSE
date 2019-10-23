# ===================================================================
# The contents of this file are dedicated to the public domain.  To
# the extent that dedication to the public domain is not available,
# everyone is granted a worldwide, perpetual, royalty-free,
# non-exclusive license to exercise all rights associated with the
# contents of this file for any purpose whatsoever.
# No rights are reserved.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
# ===================================================================

'''
Created on Fri Oct 19 21:39:56 2019


@author James Harle

$Last commit on:$
'''

import matplotlib.cm as cm
import numpy as np
import matplotlib.pyplot as plt
import NEMO_IO_utils as nemo_io
from mpl_toolkits.mplot3d.axes3d import Axes3D
import seaborn as sns
from matplotlib.colors import LightSource
sns.set(color_codes=True)

class TestConfig(object):
    """
    An object that stores the information used in constructing a idealised
    NEMO configuration (with or without an AGRIF zoom)
    """
    def __init__(self, zco=False, zps=False, sco=False, isfcav=False):
        """ 
        Initialise a NEMO test case object.

        Args:

        Returns:

        """





nmatch = 2
        # Parent
        self.zco          = zco
        self.zps          = zps
        self.sco          = sco 
        self.isfcav       = isfcav
        # AGRIF
        self.1_zco        = False
        self.1_zps        = False
        self.1_sco        = False 
        self.1_isfcav     = False
        self.zoom_ist     = None
        self.zoom_ien     = None
        self.zoom_jst     = None
        self.zoom_jen     = None
        # General
        self.bmatch       = True    #TODO: check this is ok as default
        self.nghost       = 3       #TODO: check this is ok as default

    def set_config(self, cfg):
        
        
    def set_agrif(self, iraf, jraf)
        
    def _calc_rmax(dept):
        ''' 
        R-max: measure of steepness
    
        This function returns the slope steepness criteria R-max, which is simply
        (H[0] - H[1]) / (H[0] + H[1])
    
        Args:
            zt    (np.ndarray): depth array at t-points (m)
            
        Returns:
            rm    (np.ndarray): r-max array at t-points (m)
        '''
        rmax_x = np.zeros_like(dept)
        rmax_y = np.zeros_like(dept)
        
        rmax_x[:,1:-1] = 0.5 * ( np.diff(dept[:,  :-1], axis=1) / 
                               (dept[:,  :-2] + dept[:, 1:-1]) +
                                 np.diff(dept[:, 1:  ], axis=1) / 
                               (dept[:, 1:-1] + dept[:, 2:  ]) )
        rmax_y[1:-1,:] = 0.5 * ( np.diff(dept[ :-1, :], axis=0) / 
                               (dept[ :-2, :] + dept[1:-1, :]) +
                                 np.diff(dept[1:  , :], axis=0) / 
                               (dept[1:-1, :] + dept[2:  , :]) )
        
        rmax = np.maximum(np.abs(rmax_x), np.abs(rmax_y))
        
        return rmax
    
    def _calc_stiff(dept):
        ''' 
    
        Args:
            zt    (np.ndarray): depth array at t-points (m)
            
        Returns:
            sf    (np.ndarray): stiffness factor array at t-points (m)
        '''

        sf = 
        
        return sf
    
    def channel(lx = 390., ly = 294., dx = 6., dy = 6., stiff = 1, H0 = 4500.):
        ''' 
        Creates a channel idealised domain
    
        This function produces bathymetry of a channel with a Gaussian 
        seamount in order to simulate an idealised test case.
    
        Args:
            lx, ly     (float): dimensions of the domain (km)
            dx, dy     (float): length of grid cells (km)
            stiff      (float): scale factor for steepness of seamount
    
        Returns:
            zt    (np.ndarray): depth array at t-points (m)
            rm    (np.ndarray): r-max array at t-points (m)
            x, y  (np.ndarray): coordinates (km)
        '''
        
        # Work out the number of grid points in I/J
        nx, ny = np.rint(lx/dx), np.rint(ly/dy)
        nx_int, ny_int = nx.astype(np.int32), ny.astype(np.int32)
        
        # Redefine cell size to factor into lx and ly
        dx, dy = lx/nx, ly/ny
        
        # Define T-grid
        x = np.linspace(0. + dx/2., lx - dx/2., num=nx_int)
        y = np.linspace(0. + dy/2., ly - dy/2., num=ny_int)
        x, y = np.meshgrid(x, y)
    
        # Define Gaussian bowl and R-max
        zt = np.ma.array(H0*(1. - 0.9*np.exp(-(stiff/40.**2 * ((x-lx/2.)**2 
                                      + (y-ly/2.)**2)))), mask=False)
        rm = np.ma.array(_calc_rmax(zt), mask=False)
        
        # Predfine default zoom incase AGRIF activated
        self.zoom_ist = 10
        self.zoom_ien = 40
        self.zoom_jst = 1
        self.zoom_jen = 2
        
        return x, y, zt, rm
    
    def slope(lx  = 500., ly  = 500., dx = 10., dy = 10., 
              lx2 =  50., ly2 = 100., stiff = 1):
        ''' 
        Creates a idealised sloped domain
    
        This function produces sloped bathymetry in order to 
        simulate an idealised overflow test case.
    
        Args:
            lx, ly     (float): dimensions of the domain (km)
            dx, dy     (float): length of grid cells (km)
            lx2, ly2   (float): dimensions of dense water inlet (km)
            stiff      (float): scale factor for steepness of seamount
    
        Returns:
            zt    (np.ndarray): depth array at t-points (m)
            rm    (np.ndarray): r-max array at t-points (m)
            x, y  (np.ndarray): coordinates (km)
        '''
        
        # Work out the number of grid points in I/J
        nx, ny = np.rint(lx/dx), np.rint(ly/dy)
        nx_int, ny_int = nx.astype(np.int32), ny.astype(np.int32)
        my = (ny_int+1)/2
        
        # Redefine cell size to factor into lx and ly
        dx, dy = lx/nx, ly/ny
      
        # Add the ledge to nx
        nx2, ny2 = np.floor(lx2/dx), np.floor(ly2/dy)
        nx2_int, ny2_int = nx2.astype(np.int32), ny2.astype(np.int32)
        lx2 = nx2*dx
        
        # Define T-grid
        x = np.linspace(-lx2 + dx/2., lx - dx/2., num=nx_int + nx2_int)
        y = np.linspace(   0.+ dy/2., ly - dy/2., num=ny_int)
        x, y = np.meshgrid(x, y)
    
        # Define Gaussian bowl and R-max
        zt = np.ma.array(stiff * 10. * x + 600., mask=False)
        zt = np.ma.where(zt< 600.,  600., zt)
        zt = np.ma.where(zt>3600., 3600., zt)
        rm = np.ma.array(_calc_rmax(zt), mask=False)
        
        # Create source region (i.e. channel where dense water will be released)
        zt.data[:my,:nx2_int] = 0.
        zt.data[my+ny2_int:,:nx2_int] = 0.
        
        return x, y, zt, rm
    
    def half_bowl(lx  = 500., ly  = 500., dx = 10., dy = 10., 
                  lx2 =  50., ly2 = 100., stiff = 1):
        ''' 
        Creates a half bowl idealised domain
    
        This function produces a half Gaussian-like bathymetry 
        in order to simulate an idealised overflow test case.
    
        Args:
            lx, ly     (float): dimensions of the domain (km)
            dx, dy     (float): length of grid cells (km)
            lx2, ly2   (float): dimensions of dense water inlet (km)
            stiff      (float): scale factor for steepness of slope
    
        Returns:
            zt    (np.ndarray): depth array at t-points (m)
            rm    (np.ndarray): r-max array at t-points (m)
            x, y  (np.ndarray): coordinates (km)
        '''
        
        # Work out the number of grid points in I/J
        nx, ny = np.rint(lx/dx), np.rint(ly/dy)
        nx_int, ny_int = nx.astype(np.int32), ny.astype(np.int32)
        mx, my = (nx_int+1)/2, (ny_int+1)/2
        
        # Redefine cell size to factor into lx and ly
        dx, dy = lx/nx, ly/ny
      
        # Add the ledge to nx
        nx2, ny2 = np.floor(lx2/dx), np.floor(ly2/dy)
        nx2_int, ny2_int = nx2.astype(np.int32), ny2.astype(np.int32)
        lx2 = nx2*dx
        
        # Define T-grid
        x = np.linspace(-lx2 + dx/2., lx - dx/2., num=nx_int + nx2_int)
        y = np.linspace(   0.+ dy/2., ly - dy/2., num=ny_int)
        x, y = np.meshgrid(x, y)
    
        # Define Gaussian bowl and R-max
        zt = np.ma.array(3000.*np.exp(-(stiff * 10./(lx/2.)**2)*(x-lx/2.)**2 
                                      -(stiff * 10./(ly/2.)**2)*(y-ly/2.)**2) 
                        + 600., mask=False)
        rm = np.ma.array(_calc_rmax(zt), mask=False)
        
        # Open bowl
        zt.data[:,mx+nx2_int:] = zt.data[:,mx+nx2_int-1:mx+nx2_int]
        rm.data[:,mx+nx2_int:] = rm.data[:,mx+nx2_int-1:mx+nx2_int]
        
        # Create source region (i.e. channel where dense water will be released)
        zt.data[:my,:nx2_int] = 0.
        zt.data[my+ny2_int:,:nx2_int] = 0.
        
        return x, y, zt, rm
    
    def plot_dom(x, y, zt, rm, title):
        ''' 
        Plot idealised domain
    
        Plots the idealised overflow test case as a 3D surface
        plot, coloured with the associated r-max value for reference.
    
        Args:
            x, y  (np.ndarray): coordinates (km)
            zt    (np.ndarray): depth array at t-points (m)
            rm    (np.ndarray): r-max array at t-points (m)
            
        Returns:
        '''
        
        # Set figure and axes handles
        fig = plt.figure(figsize=plt.figaspect(0.33)*1.5)
        ax  = fig.gca(projection='3d')
    
        # Create light source object, shade data of r-max and rgb array
        ls = LightSource(azdeg=0, altdeg=90)
        rgb = ls.shade(rm, plt.cm.RdYlBu)
    
        # Plot surface
        ax.plot_surface(x, y, zt, rstride=1, cstride=1, linewidth=0,
                               antialiased=False, facecolors=rgb)
    
        # Scrape colour map data for r-max
        m = cm.ScalarMappable(cmap=cm.RdYlBu)
        m.set_array(rm)
    
        # Tidy and annotate
        cb = fig.colorbar(m,shrink=0.8, aspect=8)
        ax.set_zlim(0,)
        ax.invert_zaxis()
        cb.set_label('R-max')
        ax.set_xlabel('X [km]')
        ax.set_ylabel('Y [km]')
        ax.set_zlabel('Depth [m]') 
        
        dx, dy = x[0,1]-x[0,0], y[1,0]-y[0,0]
        title = title+' (dx={}, dy={})'.format(dx, dy)
        plt.title(title,fontweight='bold')
        plt.show()
        
    def child_gen(x, y, zt, rm, title):
        ''' 
        Plot idealised domain
    
        Plots the idealised overflow test case as a 3D surface
        plot, coloured with the associated r-max value for reference.
    
        Args:
            x, y  (np.ndarray): coordinates (km)
            zt    (np.ndarray): depth array at t-points (m)
            rm    (np.ndarray): r-max array at t-points (m)
            
        Returns:
        '''
    
    def _coord_gen(self, dx, dy):
        ''' 
        Populate the coordinate arrays given a pre-defined T-Grid.
    
        Args:
            dx, dy  (float): coordinates (km)
            
        Returns:
        '''
        lonv = lont
        lonf = lonu
        latu = latt
        latf = latv
    
        self.e1t = np.ones_like(self.lont) * dx
        self.e2t = np.ones_like(self.lont) * dy
        self.e1u = np.ones_like(self.lont) * dx
        self.e2u = np.ones_like(self.lont) * dy
        self.e1v = np.ones_like(self.lont) * dx
        self.e2v = np.ones_like(self.lont) * dy
        self.e1f = np.ones_like(self.lont) * dx
        self.e2f = np.ones_like(self.lont) * dy