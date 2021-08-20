# pip install wheel
# pip install pipwin
# pipwin install numpy
# pipwin install pandas
# pipwin install shapely
# pipwin install gdal
# pipwin install fiona
# pipwin install pyproj
# pipwin install six
# pipwin install rtree
# pipwin install geopandas
# https://www.naturalearthdata.com/downloads/


#  Import libraries ----------------------------------------------------------------------------
import geopandas as gpd
import matplotlib.pyplot as plt
import matplotlib as mpl
import serial
import time
import numpy as np
from tkinter import * 
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from tkinter import colorchooser
import os
import urllib.request
from bs4 import BeautifulSoup as sp
from mpl_toolkits.axes_grid1 import make_axes_locatable
from geofeather.core import from_geofeather
from shapely.geometry.polygon import Point


# Get date from sst file name
def get_date(txt):
    date = txt.split('.')[5]
    yy = date[0:4]
    mm = date[4:6]
    dd = date[6:8]
    date = mm+'/'+dd+'/'+yy
    return date

# Web scrapper
# https://pythonprogramminglanguage.com/get-links-from-webpage/
# https://www.cpc.ncep.noaa.gov/products/GIS/GIS_DATA/sst_oiv2/
def get_sst_links():
    global link_dates, shp_links
    link_dates = []
    n=0
    try:
        myurl = "https://www.cpc.ncep.noaa.gov/products/GIS/GIS_DATA/sst_oiv2/"
        client = urllib.request.urlopen(myurl)
        html_page = client.read()
        client.close()
        sp_page = sp(html_page, "html.parser")
        links = []
        shp_links = []
        for link in sp_page.findAll('a'):
            links.append(link.get('href'))
        
        for item in links:
            txt = str(item)
            if 'ftp' and 'shp' in txt:
                date = get_date(txt)
                link_dates.append(date) # same index as shp_links
                shp_links.append(txt)
                print(date)
                print(txt)
                n = n+2
                txt = date + ":\n" + txt + "\n"
                link_position = str(n)+ '.0'
                tag_position = str(n-1)+ '.0'
                date_end = str(n-1)+ '.10'
                tag_name = 'dte'+str(n)
                links_text.insert(link_position, txt)
                links_text.tag_add(tag_name,tag_position,date_end)
                links_text.tag_config(tag_name, background='black', foreground='white', underline=True)
                
                
    except:
        links_text.insert('1.0', 'ERROR: Check your internet connection')
    
    return n/2


# Download SST file
# https://www.tutorialspoint.com/downloading-files-from-web-using-python
def download_file():
    global sst_dir1
    sst_dir1 = ''
    try:
        date_selected = dates.get()
        i = link_dates.index(date_selected) # index of matching date is equal to link index
        print( i, date_selected, shp_links[i])
        mylink = shp_links[i]
        if mylink.find('/'):
            sst_filename = mylink.rsplit('/',1)[1]
            print("filename:", sst_filename)
        # dest_dir = 'world_map_files/SST/' + sst_filename
        cd = os.getcwd()
        print(cd) # C:\Users\christian.saiz\Documents\Coding_courses\Python
        if "C:\\" in cd:
            cd = cd.replace("C:\\",'')
            print(cd)
        sst_file_dir = cd + "\\world_map_files\\SST\\" + sst_filename
        print("downloaded file:", sst_file_dir)
        urllib.request.urlretrieve(mylink, "C:\\" + sst_file_dir)
        download_status.config(background='green', foreground='white')
        txt = sst_filename + ' downloaded' 
        downloaded.set(txt)
    except:
        print("Couldn't download the file")
        download_status.config(background='red', foreground='white')
        downloaded.set('Download error')
    sst_dir1 = sst_file_dir

def browse_sst_file(): 
    global sst_dir2
    sst_dir2 = ''
    try:
        sst_file_dir = filedialog.askopenfilename(title ="Select sst shapefile zip to load...")
        if "C:/" in sst_file_dir:
            sst_file_dir = sst_file_dir.replace("C:/",'')
        fselected = sst_file_dir.split('/')[-1]
        print("SST file browsed:", sst_file_dir, fselected)
        # C:/Users/christian.saiz/Documents/Coding_courses/Python/world_map_files/SST/sst_io.20210509.shp
        download_status.config(background='green', foreground='white')
        txt = fselected + ' selected' 
        downloaded.set(txt)
        sst_dir2 = sst_file_dir
    except:
        print("Couldn't select the file")
        download_status.config(background='red', foreground='white')
        downloaded.set('Selection error')


# Sea Surface Temperature
# https://www.cpc.ncep.noaa.gov/products/GIS/GIS_DATA/sst_oiv2/index.php
def load_sst():
    try:
        dir = sst_dir1
    except:
        dir = sst_dir2
    dir.replace('\\',"/")
    print("Plotting dir:", dir)
    sst = gpd.read_file('zip:///'+ dir)
    # clip area if !0
    if(radius.get()=='0'):
        pass
    else:
        sst = gpd.clip(sst,area)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right",size ="2%", pad=0.1)
    sst.plot(ax = ax, column = 'Contour', cmap = sst_color.get(), legend = True, cax=cax )
    cax.legend(title='SST [°C]')
    cax.tick_params(labelsize=10)
    cax.tick_params(grid_markevery = 2.0)
    cax.locator_params(nbins=20)
    plt.sca(ax) # set current axe back to ax, not cax for label
# valid keywords are ['size', 'width', 'color', 'tickdir', 'pad', 'labelsize', 'labelcolor', 'zorder', 'gridOn', 'tick1On', 'tick2On', 'label1On', 'label2On', 'length', 'direction', 'left', 'bottom', 'right', 'top', 'labelleft', 'labelbottom', 'labelright', 'labeltop', 'labelrotation', 'grid_agg_filter', 'grid_alpha', 'grid_animated', 'grid_antialiased', 'grid_clip_box', 'grid_clip_on', 'grid_clip_path', 'grid_color', 'grid_contains', 'grid_dash_capstyle', 
# 'grid_dash_joinstyle', 'grid_dashes', 'grid_data', 'grid_drawstyle', 'grid_figure', 'grid_fillstyle', 'grid_gid', 'grid_in_layout', 'grid_label', 'grid_linestyle', 'grid_linewidth', 'grid_marker', 'grid_markeredgecolor', 'grid_markeredgewidth', 'grid_markerfacecolor', 'grid_markerfacecoloralt', 'grid_markersize', 'grid_markevery', 'grid_path_effects', 'grid_picker', 'grid_pickradius', 'grid_rasterized', 'grid_sketch_params', 'grid_snap', 'grid_solid_capstyle', 'grid_solid_joinstyle', 'grid_transform', 'grid_url', 'grid_visible', 'grid_xdata', 'grid_ydata', 'grid_zorder', 'grid_aa', 'grid_c', 'grid_ds', 'grid_ls', 'grid_lw', 'grid_mec', 'grid_mew', 'grid_mfc', 'grid_mfcalt', 'grid_ms']

# Select target file
def browse_file():
    global target_fdir
    global target_fname
    target_fdir = filedialog.askopenfilename(title ="Select target file to load...")
    fname_start = target_fdir.rfind("/") # return the index of the last \
    target_fname = target_fdir[target_fdir.rfind("/") + 1:] #store string from last slash to the end   
    fname.set(target_fname) 
    print('target file dir and name: ', target_fdir,target_fname)

# color scale for ocean -------------------------------------------------------

def ocol(x):
    global bcmap
    bcmap = plt.cm.get_cmap(bathy_color.get())
    return bcmap(x/10000)
#  End colormap ---------------------------------------------------------------


#  Land 
def load_land():
    
    # Read land files
    land = gpd.read_file(path1 + 'ne_10m_land.shp')
    land.crs = "EPSG:4326" 
    # land = land.set_crs("EPSG:4326")
    misland = gpd.read_file(path1 + 'ne_10m_minor_islands.shp')
    misland.crs = "EPSG:4326" 
    mislandcst = gpd.read_file(path1 + 'ne_10m_minor_islands_coastline.shp')
    mislandcst.crs = "EPSG:4326" 
    border = gpd.read_file(path1 + 'ne_10m_admin_0_boundary_lines_land.shp')
    border.crs = "EPSG:4326" 
    reefs = gpd.read_file(path1 + 'ne_10m_reefs.shp')
    reefs.crs = "EPSG:4326" 
    # clip area if !0
    if(radius.get()=='0'):
        pass
    else:
        land = gpd.clip(land,area)
        misland = gpd.clip(misland,area)
        mislandcst = gpd.clip(mislandcst,area)
        border = gpd.clip(border,area)
        reefs = gpd.clip(reefs,area)

    # Plot land features
    if (ch_sst.get() == True):
        if not land.empty:
            land.plot( ax = ax, color = 'chocolate', edgecolor = 'white', lw = 1)
        if not misland.empty:
            misland.plot( ax = ax, color = 'chocolate', edgecolor = 'white', lw = 1)
        if not mislandcst.empty:
            mislandcst.plot( ax = ax, color = 'chocolate', edgecolor = 'white', lw = 1)
        if not reefs.empty:
            reefs.plot( ax = ax, color = 'chocolate', edgecolor = 'white', lw = 0.5)
    else:
        if not land.empty:
            land.plot( ax = ax, color = 'chocolate', edgecolor = 'black', lw = 0.1)
        if not misland.empty:
            misland.plot( ax = ax, color = 'chocolate', edgecolor = 'black', lw = 0.1)
        if not mislandcst.empty:
            mislandcst.plot( ax = ax, color = 'chocolate', edgecolor = 'black', lw = 0.1)
        if not reefs.empty:
            reefs.plot( ax = ax, color = 'chocolate',edgecolor = 'black', lw = 0.1)
    if not border.empty:
        border.plot( ax = ax, color = 'white', lw = 0.5)
        

# Lakes and rivers
def load_lakes_rivers():
    if(ch_rivers.get()):
        lakes_eu = gpd.read_file(path1 + 'ne_10m_lakes_europe.shp')
        lakes_na = gpd.read_file(path1 + 'ne_10m_lakes_north_america.shp')
        lakes = gpd.read_file(path1 + 'ne_10m_lakes.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            lakes_eu = gpd.clip(lakes_eu,area)
            lakes_na = gpd.clip(lakes_na,area)
            lakes = gpd.clip(lakes,area)

        if not lakes_eu.empty:
            lakes_eu.plot( ax = ax, color = 'blue', lw = 0.5)
        if not lakes_na.empty:
            lakes_na.plot( ax = ax, color = 'blue', lw = 0.5)
        if not lakes.empty:
            lakes.plot( ax = ax, color = 'blue', lw = 0.5)


    if(ch_lakes.get()):
        rivers_eu = gpd.read_file(path1 + 'ne_10m_rivers_europe.shp')
        river_na = gpd.read_file(path1 + 'ne_10m_rivers_north_america.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            rivers_eu = gpd.clip(rivers_eu,area)
            river_na = gpd.clip(river_na,area)
        if not rivers_eu.empty:
            rivers_eu.plot( ax = ax, color = 'blue', lw = 0.5)
        if not river_na.empty:
            river_na.plot( ax = ax, color = 'blue', lw = 0.5)

# Ocean bathy LD
def load_bathy():
    # Read ocean files
    bathy0 = gpd.read_file(path1 + 'ne_10m_bathymetry_L_0.shp')
    bathy200 = gpd.read_file(path1 + 'ne_10m_bathymetry_K_200.shp')
    bathy1k = gpd.read_file(path1 + 'ne_10m_bathymetry_J_1000.shp')
    bathy2k = gpd.read_file(path1 + 'ne_10m_bathymetry_I_2000.shp')
    bathy3k = gpd.read_file(path1 + 'ne_10m_bathymetry_H_3000.shp')
    bathy4k = gpd.read_file(path1 + 'ne_10m_bathymetry_G_4000.shp')
    bathy5k = gpd.read_file(path1 + 'ne_10m_bathymetry_F_5000.shp')
    bathy6k = gpd.read_file(path1 + 'ne_10m_bathymetry_E_6000.shp')
    bathy7k = gpd.read_file(path1 + 'ne_10m_bathymetry_D_7000.shp')
    bathy8k = gpd.read_file(path1 + 'ne_10m_bathymetry_C_8000.shp')
    bathy9k = gpd.read_file(path1 + 'ne_10m_bathymetry_B_9000.shp')
    bathy10k = gpd.read_file(path1 + 'ne_10m_bathymetry_A_10000.shp')
    bathy0.crs = "EPSG:4326" 
    bathy200.crs = "EPSG:4326"
    bathy1k.crs = "EPSG:4326"
    bathy2k.crs = "EPSG:4326"
    bathy3k.crs = "EPSG:4326"
    bathy4k.crs = "EPSG:4326"
    bathy5k.crs = "EPSG:4326"
    bathy6k.crs = "EPSG:4326"
    bathy7k.crs = "EPSG:4326"
    bathy8k.crs = "EPSG:4326"
    bathy9k.crs = "EPSG:4326"
    bathy10k.crs = "EPSG:4326"
    # clip area if !0
    if(radius.get()=='0'):
        pass
    else:
        bathy0 = gpd.overlay(bathy0, area, how ='intersection')
        bathy200 = gpd.overlay(bathy200, area, how ='intersection')
        bathy1k = gpd.overlay(bathy1k, area, how ='intersection')
        bathy2k = gpd.overlay(bathy2k, area, how ='intersection')
        bathy3k = gpd.overlay(bathy3k, area, how ='intersection')
        bathy4k = gpd.overlay(bathy4k, area, how ='intersection')
        bathy5k = gpd.overlay(bathy5k, area, how ='intersection')
        bathy6k = gpd.overlay(bathy6k, area, how ='intersection')
        bathy7k = gpd.overlay(bathy7k, area, how ='intersection')
        bathy8k = gpd.overlay(bathy8k, area, how ='intersection')
        bathy9k = gpd.overlay(bathy9k, area, how ='intersection')
        bathy10k = gpd.overlay(bathy10k, area, how ='intersection')

    # Plot ocean features
    if not bathy0.empty:
        bathy0.plot( ax = ax, color = ocol(0))
    if not bathy200.empty:
        bathy200.plot( ax = ax, color = ocol(200))
    if not bathy1k.empty:
        bathy1k.plot( ax = ax, color = ocol(1000))
    if not bathy2k.empty:
        bathy2k.plot( ax = ax, color = ocol(2000))
    if not bathy3k.empty:
        bathy3k.plot( ax = ax, color = ocol(3000))
    if not bathy4k.empty:
        bathy4k.plot( ax = ax, color = ocol(4000))
    if not bathy5k.empty:
        bathy5k.plot( ax = ax, color = ocol(5000))
    if not bathy6k.empty:
        bathy6k.plot( ax = ax, color = ocol(6000))
    if not bathy7k.empty:
        bathy7k.plot( ax = ax, color = ocol(7000))
    if not bathy8k.empty:
        bathy8k.plot( ax = ax, color = ocol(8000))
    if not bathy9k.empty:
        bathy9k.plot( ax = ax, color = ocol(9000))
    if not bathy10k.empty:
        bathy10k.plot( ax = ax, color = ocol(10000))

# Ocean bathy HD
# https://opendem.info/download_bathymetry.html
def load_bathyhd():
    global bathyhd_file
    bathyhd_file = from_geofeather('world_map_files\\HDbathy_converted\\hdbathy.feather') # method to read gdf faster
    if(radius.get()=='0'):
        pass
    else:
        filter_index = bathyhd_file.within(area.loc[0,'geometry'])
        bathyhd_file = bathyhd_file.loc[filter_index]
    #bathyhd_file = bathyhd_file.cx[-80.0:-70.0, 20.0:30.0] # filter and show only shapes within the boundary box
    # bathyhd_file = gpd.read_file(path3 + 'gebco_derived_polygons.shp') # this is slower
    # bathyhd_file.plot( ax = ax, column = 'gridcode', cmap = 'ocean', legend = True ) # this worked with original file
    bathyhd_file.plot( ax = ax, column = 'Depth', cmap = bathy_color.get(),edgecolor = 'black', lw = 0.1, legend = True ) # modified shp includes Depth column


# Aux map lines
def load_grid():
    print('grid: ',ch_grid30.get(),ch_grid20.get(),ch_grid15.get(),ch_grid10.get(),ch_grid5.get(),ch_grid1.get())
    if(ch_grid30.get()):
        grid30 = gpd.read_file(path1 + 'ne_10m_graticules_30.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            grid30 = gpd.clip(grid30,area)
        grid30.plot( ax = ax, color = 'red', lw = 0.25) 

    if(ch_grid20.get()):
        grid20 = gpd.read_file(path1 + 'ne_10m_graticules_20.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            grid20 = gpd.clip(grid20,area)
        grid20.plot( ax = ax, color = 'red', lw = 0.25) 

    if(ch_grid15.get()):
        grid15 = gpd.read_file(path1 + 'ne_10m_graticules_15.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            grid15 = gpd.clip(grid15,area)
        grid15.plot( ax = ax, color = 'red', lw = 0.25) 

    if(ch_grid10.get()):
        grid10 = gpd.read_file(path1 + 'ne_10m_graticules_10.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            grid10 = gpd.clip(grid10,area)
        grid10.plot( ax = ax, color = 'red', lw = 0.25)

    if(ch_grid5.get()):
        grid5 = gpd.read_file(path1 + 'ne_10m_graticules_5.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            grid5 = gpd.clip(grid5,area)
        grid5.plot( ax = ax, color = 'red', lw = 0.25)

    if(ch_grid1.get()):
        grid1 = gpd.read_file(path1 + 'ne_10m_graticules_1.shp')
        # clip area if !0
        if(radius.get()=='0'):
            pass
        else:
            grid1 = gpd.clip(grid1,area)
        grid1.plot( ax = ax, color = 'red', lw = 0.25)
    


# Spherical Law of Cosines - distance between 2 locations
def dist2target(lat1, lon1, lat2, lon2):
    p1 = np.deg2rad(lat1) # lat1
    p2 = np.deg2rad(lat2) # lat2
    dg = np.deg2rad( lon2 - lon1 ) # delta long
    R = 6371E3 # Earth radius in m
    d = np.arccos(np.sin(p1)*np.sin(p2)+np.cos(p1)*np.cos(p2)*np.cos(dg))*R
    return int(d/1000) # in km

# Draw clicked points
def draw_point(xd,yd, fsize, fcolor, msize, text, mark ):
    plt.plot(xd,yd, mark, color = fcolor, markersize = msize)
    plt.text(xd,yd+0.01, text , color = fcolor, fontsize = fsize, weight = 'bold')
# Endof draw_location -----------------------------------------------------------------------------------

# Country names ------------------------------------------------------------------------------------------
def load_names(fname):
    cnfile = open (path1 + fname, 'r')
    country_list = cnfile.readlines()
    for line in country_list:
        if 'END' in line:
            break
        try:
            line_cols = line.split()
            # print(line_cols)
            gdf_point = gpd.GeoDataFrame( [ {'geometry': Point(float(line_cols[1]),float(line_cols[0])) } ] )
            gdf_point.crs = "EPSG:4326"
            # print( gdf_point.within(area)[0] )
            
            if( gdf_point.within(area)[0]  == True):
                print(line_cols[2], 'is included')
                plt.text(float(line_cols[1]),float(line_cols[0]), line_cols[2] , color = 'black', fontsize = 6)
        except:
            pass
            
    cnfile.close()
# End load_names ----------------------------------------------------------------------------------------

# Targets file ------------------------------------------------------------------------------------------
# Text file format is LAT LON NAME COLOR FONT_SIZE
# Use # for comments in the text file and end the file filling END
def load_targets():   
    try: 
        tgtfile = open(target_fdir, 'r')
    except:
        tgtfile = open('targets.txt', 'r')

    target_list = tgtfile.readlines()
    for t_line in target_list:
        if 'END' in t_line:
            break
        if '#' in t_line:
            continue
        try:
            t_line_cols = t_line.split()
            dist = dist2target(float(t_line_cols[0]), float(t_line_cols[1]), gps_val[0], gps_val[1])
            print(t_line, dist)
            if (dist > 1.0):
                target_text = t_line_cols[2] + "\n" + "{:.3f}".format(float(t_line_cols[0])) + ' ' + "{:.3f}".format(float(t_line_cols[1])) + '\nDist: ' + str(dist) + 'km'
            else:
                target_text = t_line_cols[2] + "\n" + "{:.3f}".format(float(t_line_cols[0])) + ' ' + "{:.3f}".format(float(t_line_cols[1])) + '\nDist: ' + str(dist*1000) + 'm'
            
            draw_point(float(t_line_cols[1]), float(t_line_cols[0]), t_line_cols[4], t_line_cols[3], 3, target_text, 'o')
        except:
            pass
    tgtfile.close()
#  End load_targets ------------------------------------------------------------------------------------

# Points of interest -------------------------------------------------------
# https://matplotlib.org/stable/users/event_handling.html
def onclick(event):
    xf,yf = event.x, event.y
    xd,yd = event.xdata,event.ydata
    btn = event.button
    print(btn, 'xy_figure:',xf,yf,' - ', 'xy_data:', xd,yd)
    dist = dist2target(yd, xd, gps_val[0], gps_val[1])
    if(dist > 1.0):
        point_text = 'Lat: ' + "{:.3f}".format(yd) + ' Lon: ' + "{:.3f}".format(xd) + '\nDist: ' + str(dist) + 'km'
    else:
        point_text = 'Lat: ' + "{:.3f}".format(yd) + ' Lon: ' + "{:.3f}".format(xd) + '\nDist: ' + str(dist*1000) + 'm'

    if(btn == 2): # 1 LEFT, 2 MIDDLE, 3 RIGHT, 8 wheel back, 9 wheel forward
        draw_point(xd,yd, 8, point_color[1], 3, point_text, 'X') # position in map, fontsize, color, marker size, text, mark
        print('margins: ', ax.margins())
        plt.show()   

# End points of interest ----------------------------------------------------- 


def color_pick():
    global point_color
    point_color = colorchooser.askcolor(title = "Choose color for points of interest")
    print(point_color)
    color_lbl.configure(background=point_color[1])

# GPS ---------------------------------------------------------------------------------------------------
class gpstime:
    def __init__(self, hh, mm, ss):
        self.hh = hh
        self.mm = mm
        self.ss = ss

class gpsdate:
    def __init__(self, d, m, y):
        self.d = d
        self.m = m
        self.y = y

def deg(x):
    lx = len(str(int(float(x))))
    xdeg = str(int(float(x)))
    if ( lx == 3):
        return xdeg[0:1]
    elif ( lx == 4):
        return xdeg[0:2]
    elif ( lx == 5):
        return xdeg[0:3]

def minu(x):
    lx = len(str(int(float(x))))
    xmin = str(float(x))
    if (lx == 3):
        return xmin[1:]
    elif (lx == 4):
        return xmin[2:]
    elif (lx == 5):
        return xmin[3:]    

def lat_deg2dec(lat_deg,lat_min, ns):
    a = 1
    if ns == 'N': a = 1
    if ns == "S": a = -1
    return (lat_deg + lat_min/60)*a

def lon_deg2dec( lon_deg,lon_min, ew):
    a = 1
    if ew == 'E': a = 1
    if ew == "W": a = -1
    return (lon_deg + lon_min/60)*a


# GPS serial port
def load_gps():
    global gps_val, gps_text
    try:
        dt = 2
        print('Loading GPS port/baud: ', sport.get(), sbaud.get())
        myport = serial.Serial(port = sport.get(), baudrate = sbaud.get(), bytesize=8,timeout=2, stopbits=serial.STOPBITS_ONE)
        t0 = time.time()
        while( (time.time() - t0) < dt):
                if(myport.in_waiting > 0):
                    line = myport.readline()
                    dline = line.decode('Ascii')
                    if( '$GPRMC' in dline):
                        lline = dline.split(',')
                        utc_date = lline[9]
                        utc_time = lline[1]
                        utc_date = gpsdate(utc_date[0:2], utc_date[2:4], utc_date[4:6])
                        utc_time = gpstime(utc_time[0:2], utc_time[2:4], utc_time[4:6])
                        lat = lline[3]
                        lat_deg = int(deg(lat))
                        lat_min = float(minu(lat))
                        ns = lline[4]
                        lon = lline[5]
                        lon_deg = int(deg(lon))
                        lon_min = float(minu(lon))
                        ew = lline[6]
                        v = lline[7]
                        hdg = lline[8]
                        txt = "Date: {}/{}/{}  Time: {}:{}:{} Lat: {}° {}' {}  Lon: {}° {}' {} Speed: {} kt  Heading: {} °"
                        gps_text = txt.format( utc_date.m, utc_date.d, utc_date.y, utc_time.hh, utc_time.mm, utc_time.ss , lat_deg, lat_min, ns, lon_deg, lon_min, ew, v, hdg )
                        print(gps_text)
                        lat_dec = lat_deg2dec(lat_deg,lat_min, ns) 
                        lon_dec = lon_deg2dec( lon_deg,lon_min, ew)
                        mypos = [float(lat_dec), float(lon_dec), float(v), float(hdg)]
                        print(mypos)
                        # gps_window = Toplevel(root)
                        # gps_window.title('GPS data')
                        

        # ttk.Label(gps_window, text = gps_text, justify = LEFT).grid(row = 0, column = 0, sticky = 'nesw')
        title_font = font = {'family': 'serif', 'color':  'green', 'weight': 'normal', 'size': 8, }
        ax.set_title('WORLD MAP: '+ gps_text, fontdict = title_font )                
        myport.close()
    except:      
        mypos = [float(dlat.get()),float(dlon.get()),float(dspeed.get()),float(dheading.get())]  
        print('Ploting position ', mypos[0],mypos[1],mypos[2],mypos[3]) 
        messagebox.showerror(title = 'GPS error', message = "Can't connect with GPS\nPlotting manual GPS values" )
     
    gps_val = mypos
     
# End GPS ------------------------------------------------------------------
   

# Color bar --------------------------------------------
def color_bar():
    bar, axbar = plt.subplots(num = 2, figsize = (2,3), squeeze=True)
    # bar.set_figheight(5)
    # axbar.set_title('OCEAN DEPTH [m]')
    bounds = [0,200,1000,2000,3000,4000,5000,6000,7000,8000,9000,10000]
    norm = mpl.colors.BoundaryNorm(bounds, bcmap.N)
    bar.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap = bathy_color.get()),
    cax=axbar, orientation='vertical', label='Ocean depth [m]', aspect = 50, spacing = 'proportional', ticks = bounds)
    plt.tight_layout()
    # https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.colorbar.html

def plot_gps():
    lat = str("{:.3f}".format(gps_val[0]))
    lon = str("{:.3f}".format(gps_val[1]))
    text = 'Lat: ' + lat + ' Lon: ' + lon + "\nv: " + str(gps_val[2]) + ' hdg: ' + str(gps_val[3])
    draw_point(gps_val[1], gps_val[0]+0.01, 8, 'magenta', 3, text, 'D')
    plt.arrow(gps_val[1],gps_val[0]+0.01, gps_val[2]*np.sin(np.deg2rad(gps_val[3]))/30, gps_val[2]*np.cos(np.deg2rad(gps_val[3]))/30 , color = 'magenta',head_width = 0.01)

def create_map_window():
    global ax, myMap
    myMap, ax = plt.subplots(num = 1, figsize = (14,7))
    myMap.tight_layout(pad = 1.2)
    myMap.subplots_adjust(hspace=0.5, wspace=0.5)
    ax.margins(0)
    # plt.subplots_adjust(left=1, right=1, top=1, bottom=1)
    ax.set_title('WORLD MAP')
    ax.set_ylabel('Latitude')
    ax.set_xlabel('Longitude')


def start_map():
    print("Button pressed")
    print('Target check is ',ch_targets.get())
    load_gps() # retrieve gps data and write to global variable gps_val = [lat,lon,v,hdg]
    extent_obj() # create circunference to plot, centered at gps location
    
    if(ch_sst.get()):
        load_sst()

    plot_gps() # plot position and arrow
    load_land() # land features
    load_lakes_rivers()
    load_grid() # 30 20 15 10 5 1 degree grids
    
    if(ch_bathy.get()):
        load_bathy() # ocean features
    if(ch_bathyhd.get()):
        load_bathyhd()
    if(ch_cnames.get()):
        load_names('Country.txt') # name of the countries
    if(ch_targets.get()):
        load_targets() # read from text file
    if(ch_bathy.get()):
        color_bar()
    plt.show()

def target_help():
    text = 'Any text file with .txt extension can be opened \nto plot any number of predefined targets on the map'
    text = text + '\n\nThe structure should be:\nLAT LON NAME COLOR FONT_SIZE\nEx. 27.255 -80.052 CTD1 blue 6\n\nAdd END at the end of the file'
    text = text + '\nUse a starting # for comments'
    messagebox.showinfo(title = 'Info about targets file...', message = text)

def gps_help():
    text = 'If a GPS serial with connection is available, configure the port\nOtherwise enter the GPS values manually'
    text = text + '\n- Lat and lon in decimal degrees\n- Speed over ground in knots\n- Heading in degrees from True North = 0°'
    text = text + '\n\nSerial parameters: bytesize = 8 / timeout = 2 / stopbits = 1'
    messagebox.showinfo(title = 'Info about gps use...', message = text)

def layers_help():
    text = 'Select the layers you want to display on the map\nPlot targets file only if available'
    text = text + '\n- Plot bathymetry layer from 0 to 10,000 m of depth\n- Plot parallels and meridians grid\n- Plot the name of all the countries on the map'
    text = text + '\n- Plot main rivers of North America and Europe\n- Plot main lakes'
    text = text + '\nUse Bathymetry LD for a faster loading and navigation'
    text = text +'\nBathymetry HD is slower, therefore a map extent of ~ 45° is recommended'
    messagebox.showinfo(title = 'Info about layers...', message = text)

def sst_help():
    text = 'The app looks up for link at the National Weather Service webpage for updated links'
    text = text + '\ncontaining new SST data, usually updated every 2 weeks.'
    text = text + '\nSelect a date from the list containing an existing sst zip file'
    text = text + '\nand then click on Select date from NWS. This automatically downloads the file'
    text = text + '\nto the root directory where the executable is located'
    text = text + '\nIf you have another sst ZIP on your folders, just browse it'
    text = text + '\nThe selected file will be shown on the status bar below the text field'
    text = text + '\nRemember to check the SST checkbox to display the SST data on the map'
    messagebox.showinfo(title = 'Info about Sea surface temperature...', message = text)

def poi_help():
    text = 'POI are points of interest to be drawn on the map by clicking anywhere'
    text = text + '\nClick on Select to choose the color\nOn the map, use the middle button of the mouse to draw a point'
    text = text + '\nThe POI will display the position and distance to the current GPs position'
    messagebox.showinfo(title = 'Info about Point of Interest (POI)...', message = text)

def bathy_help():
    text = 'There are 2 bathymetry files to plot, LD will be faster to process but less detailed than HD'
    text = text + '\n- Low Definition: resolution of 1000 m of depth'
    text = text + '\n- High Definition: resolutions ranging from: \n- 25m on shallow waters\n- 500 m on deeper waters'
    messagebox.showinfo(title = 'Info about Bathymetry...', message = text)

def colormap_help():
    text = 'Select the colormap from the following list to display the bathymetry features:'
    text = text + "\n'viridis', 'plasma', 'inferno', 'magma', 'cividis'"
    text = text + "\n'Greys', 'Purples', 'Blues', 'Greens', 'Oranges'" 
    text = text + "\n'Reds','YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu'"
    text = text + "\n'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'"
    text = text + "\n'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone', 'pink'"
    text = text + "\n'spring', 'summer', 'autumn', 'winter', 'cool', 'Wistia'"
    text = text + "\n'hot', 'afmhot', 'gist_heat', 'copper'"
    text = text + "\n'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern'"
    text = text + "\n'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg'"
    text = text + "\n\nReverse the color by adding '_r' at the end of the color name i.e. 'ocean_r'"
    messagebox.showinfo(title = 'Info about colormaps..', message = text)

def restart_all():
    plt.close('all')
    root.destroy()
    # print(os.getcwd())
    this_file = os.path.basename(__file__)
    line = 'python ' + this_file
    os.system(line)
    

def exit_all():
    plt.close('all')
    root.destroy()

# Disable sst and bathy HD selection
def cancel_sstld():
    if(ch_bathy.get() == True):
        ch_sst.set(False)
        ch_bathyhd.set(False)
        sst_btn.config(state = 'disabled')
        links_text.config(state = 'disabled')
        sst_pick.config(state = 'disabled')
        download_button.config(state = 'disabled')    
        browse_sst_button.config(state = 'disabled')  
        bathyhd_btn.config(state = 'disabled')
    else:
        sst_btn.config(state = 'normal')
        links_text.config(state = 'normal')
        sst_pick.config(state = 'normal')
        download_button.config(state = 'normal')
        browse_sst_button.config(state = 'normal')
        bathyhd_btn.config(state = 'normal')

# Disable sst and bathy LD selection
def cancel_ssthd():
    if(ch_bathyhd.get() == True):
        ch_sst.set(False)
        ch_bathy.set(False)
        sst_btn.config(state = 'disabled')
        links_text.config(state = 'disabled')
        sst_pick.config(state = 'disabled')
        download_button.config(state = 'disabled')    
        browse_sst_button.config(state = 'disabled')  
        bathy_btn.config(state = 'disabled')
    else:
        sst_btn.config(state = 'normal')
        links_text.config(state = 'normal')
        sst_pick.config(state = 'normal')
        download_button.config(state = 'normal')
        browse_sst_button.config(state = 'normal')
        bathy_btn.config(state = 'normal')

# Disable bathymetry selection
def cancel_bathy():
    if(ch_sst.get() == True):
        ch_bathy.set(False)
        ch_bathyhd.set(False)
        bathy_btn.config(state = 'disabled')
        bathyhd_btn.config(state = 'disabled')
    else:
        bathy_btn.config(state = 'normal')
        bathyhd_btn.config(state = 'normal')

def extent_obj():
    global area
    area = gpd.GeoDataFrame( [ {'geometry': Point(float(dlon.get()),float(dlat.get())).buffer(float(radius.get())) } ] )
    area.crs = "EPSG:4326"


def main():
    print('current file: ', os.path.basename(__file__))
    # global variables
    global sport, sbaud
    global dlat, dlon, dspeed, dheading
    global ch_targets, ch_bathy, ch_bathyhd, ch_cnames, ch_lakes, ch_rivers, ch_grid30, ch_grid20, ch_grid15, ch_grid10, ch_grid5, ch_grid1, ch_sst
    global fname
    global path1, path2, path3
    global root, gps_window
    global links_text
    global sst_btn, bathy_btn, bathyhd_btn, sst_pick, download_status, download_button, browse_sst_button
    global dates, downloaded, color_lbl
    global radius
    global bathy_color, sst_color

    #  Path files ----------------------------------------------------------------------------
    path0 = 'world_map_files\\'
    path1 = 'world_map_files\\10m\\'
    path2 = 'world_map_files\\SST\\'
    path3 = 'world_map_files\\HD_bathy\\'
    target_fname = 'targets.txt'
    
    # Tkinter interface
    root = Tk()
    root.title('Map settings')
    main_font = ('Courier',10,'bold')
    frame_font = ('Courier',8)
    frame0 = LabelFrame(root, text = 'Description', font = frame_font, fg = 'blue')
    frame0.grid(row = 0, column = 0, sticky = 'nesw')
    frame0.columnconfigure(1, weight = 1)
    # Image
    pic = PhotoImage(file = path0 + "aoml_logo.gif")
    small_pic = pic.subsample(15,15)
    logo = ttk.Label(frame0)
    logo.grid(row=0, column=0)   
    logo.config(image=small_pic, compound = CENTER) 
    # Description
    ttk.Label(frame0,text = "Set parameters to \nplot the World Map\n\nThen click on \nStart button",justify = CENTER, 
    anchor = CENTER, font = main_font).grid(row = 0, column = 1, sticky='nsew')
    color_lbl = ttk.Label(frame0, text = 'Select\nPOI\ncolor')
    color_lbl.grid(row = 0, column = 2, sticky='nsew')
    
    # Color picker
    color_btn = ttk.Button(frame0, text = 'Select', command = color_pick)
    color_btn.grid(row = 0, column = 3, sticky='nsew')

    # Browse target file
    frame1 = LabelFrame(root, text = 'Targets', font = frame_font, fg = 'blue')
    frame1.grid(row = 1, column = 0, sticky = 'nesw')
    frame1.columnconfigure(1, weight = 1)
    
    fname = StringVar()
    fname.set(target_fname) 
    src_label = ttk.Label(frame1, text = "Select targets file: ", justify = LEFT)
    fname_label = ttk.Label(frame1, textvariable = fname, justify = CENTER, font = ('Courier', 8, 'italic'))
    src_button = ttk.Button(frame1, text = "Browse", command = browse_file)
    src_label.grid(row = 1, column = 0,  sticky = 'w')
    fname_label.grid(row = 1, column = 1)
    src_button.grid(row = 1, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    
    # GPS configuration
    frame2 = LabelFrame(root, text = 'GPS', font = frame_font, fg = 'blue')
    frame2.grid(row = 2, column = 0, sticky = 'nesw')
    frame2.columnconfigure(1, weight = 1)
    # Enter GPS port
    ttk.Label(frame2, text = "Serial port: ", justify = LEFT).grid(row = 2, column = 0, sticky = 'w')
    sport = StringVar()
    sport.set("COM14")
    sport_entry = ttk.Entry(frame2, textvariable = sport, width = 12, justify = CENTER ).grid( row = 2, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    # Enter GPS baudrate
    ttk.Label(frame2, text = "Baud rate: ", justify = LEFT).grid(row = 3, column = 0, sticky = 'w')
    sbaud = IntVar()
    sbaud.set(4800)
    sbaud_entry = ttk.Entry(frame2, textvariable = sbaud, width = 12, justify = CENTER ).grid( row = 3, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    
    # Enter GPS default values
    dlat = StringVar()
    dlon = StringVar()
    dspeed = StringVar()
    dheading = StringVar()
    dlat.set("27.5")
    dlon.set("-80.0")
    dspeed.set("10.0")
    dheading.set("120.0")
    ttk.Label(frame2, text = "Enter values manually if GPS is not connected ", justify = LEFT, font = ('Courier', 8,'bold italic')).grid(row = 4, column = 0, columnspan = 4, sticky = 'w')
    ttk.Label(frame2, text = "Latitude dd.ddd: ", justify = LEFT).grid(row = 5, column = 0, sticky = 'w')
    lat_entry = ttk.Entry(frame2, textvariable = dlat, width = 12, justify = CENTER ).grid( row = 5, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Label(frame2, text = "Longitude dd.ddd: ", justify = LEFT).grid(row = 6, column = 0, sticky = 'w')
    lon_entry = ttk.Entry(frame2, textvariable = dlon, width = 12, justify = CENTER ).grid( row = 6, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Label(frame2, text = "Speed kt: ", justify = LEFT).grid(row = 7, column = 0, sticky = 'w')
    speed_entry = ttk.Entry(frame2, textvariable = dspeed, width = 12, justify = CENTER ).grid( row = 7, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Label(frame2, text = "Heading °: ", justify = LEFT).grid(row = 8, column = 0, sticky = 'w')
    heading_entry = ttk.Entry(frame2, textvariable = dheading, width = 12, justify = CENTER ).grid( row = 8, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    
    # Select layers
    frame3 = LabelFrame(root, text = 'Select basic layers', font = frame_font, fg = 'blue')
    frame3.grid(row = 9, column = 0, sticky = 'nesw')
    frame3.columnconfigure(1, weight = 1)
    
    ch_targets = BooleanVar()
    ch_bathy = BooleanVar()
    ch_grid30 = BooleanVar()
    ch_grid20 = BooleanVar()
    ch_grid15 = BooleanVar()
    ch_grid10 = BooleanVar()
    ch_grid5 = BooleanVar()
    ch_grid1 = BooleanVar()
    ch_cnames = BooleanVar()
    ch_lakes = BooleanVar()
    ch_rivers = BooleanVar()
    ch_bathyhd = BooleanVar()
    ttk.Checkbutton(frame3, text = 'Targets', variable = ch_targets).grid(row = 9, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    bathy_btn = ttk.Checkbutton(frame3, text = 'Bathymetry LD', variable = ch_bathy, command = cancel_sstld)
    bathy_btn.grid(row = 10, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    bathyhd_btn = ttk.Checkbutton(frame3, text = 'Bathymetry HD', variable = ch_bathyhd, command = cancel_ssthd)
    bathyhd_btn.grid(row = 10, column = 1, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Label(frame3, text = "Colormap:").grid(row = 10, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    bathy_color =  StringVar()
    bathy_color.set("ocean_r")
    ttk.Entry(frame3, textvariable = bathy_color, width = 8, justify = CENTER ).grid( row = 10, column = 3, sticky = 'nesw', padx = 2, pady = 1)

    ttk.Checkbutton(frame3, text = 'Grid 30°', variable = ch_grid30).grid(row = 11, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'Grid 20°', variable = ch_grid20).grid(row = 11, column = 1, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'Grid 15°', variable = ch_grid15).grid(row = 11, column = 2, sticky = 'w', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'Grid 10°', variable = ch_grid10).grid(row = 12, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'Grid 5°', variable = ch_grid5).grid(row = 12, column = 1, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'Grid 1°', variable = ch_grid1).grid(row = 12, column = 2, sticky = 'w', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'Country names', variable = ch_cnames).grid(row = 13, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'NA/EU Lakes', variable = ch_lakes).grid(row = 14, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Checkbutton(frame3, text = 'NA/EU Rivers', variable = ch_rivers).grid(row = 14, column = 1, sticky = 'nesw', padx = 2, pady = 1)  
    
    # SST layer configuration
    frame4 = LabelFrame(root, text = 'Sea Surface Temperature layer', font = frame_font, fg = 'blue')
    frame4.grid(row = 15, column = 0, sticky = 'nesw')
    frame4.columnconfigure(1, weight = 1)
    
    ch_sst = BooleanVar()
    sst_btn = ttk.Checkbutton(frame4, text = 'SST', variable = ch_sst, command = cancel_bathy)
    sst_btn.grid(row = 16, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    
    ttk.Label(frame4, text = "Available SST files @NOAA NWS").grid(row = 17, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    
    links_text = Text(frame4, width = 70, height = 2, font=('Courier',9))
    links_text.grid(row = 18, column = 0, columnspan = 3, sticky = 'nesw', padx = 2, pady = 1)
    text_scroll = ttk.Scrollbar(frame4, orient = VERTICAL, command = links_text.yview )
    text_scroll.grid(row = 18, column = 3, sticky = 'nse')
    links_text.config(yscrollcommand = text_scroll.set)
    
    n_links = int(get_sst_links())
    lbl = str(n_links) + ' links found'
    ttk.Label(frame4, text= lbl, font=('Courier',10, 'italic bold') ).grid(row = 17, column = 1, sticky = 'nesw', padx = 2, pady = 1)
    
    # ttk.Label(frame4, text = "Select SST date:").grid(row = 16, column = 1, sticky = 'e', padx = 2, pady = 1)
    
    dates = StringVar()
    sst_pick = Spinbox(frame4, textvariable= dates, justify=CENTER)
    sst_pick.grid(row = 16, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    sst_pick.config(values = link_dates )
    download_button = ttk.Button(frame4, text = "Select date from NWS", command = lambda: download_file())
    download_button.grid(row = 16, column = 1, sticky = 'nesw', padx = 2, pady = 1)
    browse_sst_button = ttk.Button(frame4, text = "or Browse SST zip", command = browse_sst_file)
    browse_sst_button.grid(row = 17, column = 2, sticky = 'nesw', padx = 2, pady = 1)
    downloaded = StringVar()
    downloaded.set('')
    download_status = ttk.Label(frame4, textvariable = downloaded, justify=CENTER)
    download_status.grid(row = 19, column = 0, columnspan=2, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Label(frame4, text = 'colormap:', justify=CENTER).grid(row = 19, column = 2, sticky = 'nes', padx = 2, pady = 1)
    sst_color = StringVar()
    sst_color.set('inferno')
    ttk.Entry(frame4, textvariable = sst_color, width = 8, justify = CENTER ).grid( row = 19, column = 3, sticky = 'nesw', padx = 2, pady = 1)
    
    # Map Extent
    frame5 = LabelFrame(root, text = 'Map extent', font = frame_font, fg = 'blue')
    frame5.grid(row = 20, column = 0, sticky = 'nesw')
    frame5.columnconfigure(1, weight = 1)
    ttk.Label(frame5, text = "Plot area of interest a distance R from GPS position (0° plots all)-->").grid(row = 20, column = 0, sticky = 'nesw', padx = 2, pady = 1)
    ttk.Label(frame5, text = "R [degrees dd.ddd]:").grid(row = 20, column = 1, sticky = 'nes', padx = 2, pady = 1)
    radius =  StringVar()
    radius.set("0")
    radius_entry = ttk.Entry(frame5, textvariable = radius, width = 5, justify = CENTER ).grid( row = 20, column = 2, sticky = 'nesw', padx = 2, pady = 1)

    # Start map
    start_button = ttk.Button(root, text = "Start", command = lambda: start_map()).grid(row = 21, column = 0, sticky = 'nesw', padx = 2, pady = 1) 
    
    # Create menubar
    menubar = Menu(root)
    root.config(menu = menubar)
    # Create menu item
    help_ = Menu(menubar)
    menubar.add_cascade(menu = help_ , label = 'Help')
    # Create subitems
    help_.add_command(label = 'Targets', command = target_help)
    help_.add_command(label = 'GPS', command = gps_help)
    help_.add_command(label = 'Layers', command = layers_help)
    help_.add_command(label = 'Bathymetry', command = bathy_help)
    help_.add_command(label = 'SST', command = sst_help)
    help_.add_command(label = 'POI', command = poi_help)
    
    end_map = Menu(menubar)
    menubar.add_cascade(menu = end_map, label = "Exit")
    end_map.add_command(label = 'Exit All', command = exit_all)
    end_map.add_command(label = 'Restart All', command = restart_all)
    
    create_map_window()
    add_point = myMap.canvas.mpl_connect('button_press_event',onclick)
    
    

    root.mainloop()

# Bottom code
if __name__ == "__main__":
    main()