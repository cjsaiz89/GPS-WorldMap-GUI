# GPS-WorldMap-GUI

**Description**

Global Map with Bathymetry and Sea Surface Temperature layers, GPS position and targets.

This app generates an interactive world map (or a smaller specified section), and is capable of plotting targets read from a targets file, read GPS data if connected through the serial port, select different layers such as low and high definition bathymetry, and sea surface temperature (you can either choose one from the NWS that the program downloads by itself, or select any other zip file containing sst data in .shp format (make sure it's inside a zip file).

**Files**

- *WorldMap_rev2.zip* contains
  - *WorldMap_rev2.exe* executable GUI (for Windows only)
  - *targets.txt* file example.
  - */world_map_file* contains source data to plot, such as Bathymetry HD, Bathymetry LD, borders and some SST files.

- *WorldMap_rev2.py* is the python script to run if python is installed on your computer. All libraries must be downloaded to run the script.
  - */world_map_files* directory should be located in the same directory as the python file. Extract this folder from the ZIP file.


**Map settings**
- Here you can select the color for the points of interest (POI) to draw on the map by clicking with the center button of the mouse.
Set the GPS port configuration. If there's no GPS connected but you know your position, you can just type the GPS values into the fields.

- Select the layers you want to plot:
  - *Targets* from the selected *target* file (make sure you browse and select one if you check the Targets checkbutton).
  - *Bathymetry LD* OR *Bathymetry HD* (HD will take longer to load, for an extent of R = 45 degrees, it could take around 1 min). Slect the *Colormap* you want (check *Help* menu to see different color options).
  - *Grid* lines every 30, 20, 15, 10, 5 and/or 1 degrees.
  - *Country names*.
  - *North Atlantic and Europe lakes* and/or *rivers*.
  - *SST* Sea surface Temperature.

The *SST* section automatically reads the links available at the *National Weather Service* webpage, and displays the info . If you want to download one of these, select the date from the spinbox and click on *Select date from NWS*. 
The status bar at the bottom will show that the file was *downloaded* or *selected* if browsed.
*SST* *colormap* can be selected as well as bathymetry.

**Note:** *Bathymetry* and *SST* can't be displayed simultaneously.

*Map_settings example*


![image](https://user-images.githubusercontent.com/89260258/130278901-12902520-b2ef-4ee0-948f-764c0444a427.png)

![image](https://user-images.githubusercontent.com/89260258/130278997-34629083-b466-4f6b-9dcd-184a567830f9.png)


The points of interest or *POI*, are points to draw on the map anywhere you click on with the center button of the mouse. The position and distance to the GPS position is displayed. Select the POI the color by clicking on the top right button *Select*.

*Example of POI color selection*

![image](https://user-images.githubusercontent.com/89260258/130279031-6be7421c-8b39-4925-9657-3f8bd6ced401.png)


At the bottom of the *Map settings* window, you can choose the extent of the map to plot. 
A *0* plots the entire world and for HD bathymetry in particular may take a long time to load.
Limit the area to plot around the GPS position, setting a different radius *R*.

*Example with R=30.0 degrees centered at the GPS position - SST layer*

![image](https://user-images.githubusercontent.com/89260258/130279082-0cd6b54d-f50e-4180-8ffd-66af64a50d19.png)


If GPS data is available from the serial port, the data is written as the map's title.
All maps can be zoomed in/out, panned and saved as image using the bottom left options from the plot window.

*Example of Bathymetry LD with GPS (magenta) and points of interest (red) clicked on the map*

![image](https://user-images.githubusercontent.com/89260258/130279154-73d868f3-aec2-4ede-9d52-1b56972e8ce6.png)

*Examples of Bathymetry HD*

![image](https://user-images.githubusercontent.com/89260258/130279184-c112c8de-d5ab-4ad3-8df9-dc993d60e727.png)

![image](https://user-images.githubusercontent.com/89260258/130279197-bf86d419-a594-4962-a158-32be3f1945e8.png)

![image](https://user-images.githubusercontent.com/89260258/130279203-4632a147-02a0-44fa-a73f-06d2eaae55c5.png)

![image](https://user-images.githubusercontent.com/89260258/130279210-971ae7e1-dede-4463-a210-141825633bc7.png)


**Targets file**

Create a text file with any name containing the following information:

    Lat Long Name Color Font Size
      
Rows with a *#* are ignored and can be used to comment the file.
Lat and Lon in decimal degrees (dd.ddd)


*Example of target file*

![image](https://user-images.githubusercontent.com/89260258/130279261-28b8b269-2b52-4410-b31c-89e2d660cec6.png)

*Example of target file plot*

![image](https://user-images.githubusercontent.com/89260258/130279275-f45af4c0-f44f-426f-9200-d018f062d8e4.png)

