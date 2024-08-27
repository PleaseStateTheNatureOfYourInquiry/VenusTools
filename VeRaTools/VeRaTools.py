# Author: Maarten Roos-Serote
# ORCID author: 0000 0001 5001 1347

# Version: v20240604


# Standard imports.
import numpy as np

# Custom imports.
from HandyTools import HandyTools
from DataTools import DataTools

venusSurfaceRadius = 6051.8 #km


# This is a Python class to wrangle Venus Express VeRa data.
class VeRaTools:
    '''
    This is a Python class to wrangle Venus Express VeRa data.
    '''
    

    # Read the content of a VeRa .TAB file and return a numpy array.
    @staticmethod
    def readVeRaTAB (VeRaTABFileName):
        '''
        :param VeRaTABFileName: File name (and path) of the VeRa TAB file.
        :type VeRaTABFileName: str

        :return: Values of 7 variables as a function of sounded level, numberOfLevels
        :rtype: NumPy array (7, numberOfLevels), int
        
        **Description:**
        Read the content of a VeRa .TAB file and return as a numpy array, where each row represents the values of a variable as a function of sounding level:
        
            | VeRaProfileOriginal [0] radius (km)
            | VeRaProfileOriginal [1] Temperature (K)
            | VeRaProfileOriginal [2] 1-sigma Temperature uncertainty (K)
            | VeRaProfileOriginal [3] pressure (bar)
            | VeRaProfileOriginal [4] 1-sigma pressure uncertainty (bar)
            | VeRaProfileOriginal [5] latitude (˚)
            | VeRaProfileOriginal [6] longitude (˚)
                
        '''
            
        # Read the file content as text.
        fileContentTAB = HandyTools.getTextFileContent (VeRaTABFileName)

        # Determine the number of sounded levels in the file.
        numberOfLevels = len (fileContentTAB)

        # Pascal to bar conversion factor.        
        PascalToBar = 1 / 100000

        # Extract the seven variables as a function of sounded level.        
        VeRaProfileOriginal = np.zeros ( (7, numberOfLevels) )
        for iFileLine, fileLine in enumerate (fileContentTAB):
        
            fileLineElements = [ element  for element in fileLine.split (' ')  if element != '' ]
                        
            iLevel = numberOfLevels - iFileLine - 1
            VeRaProfileOriginal [0][iLevel] = float ( fileLineElements [3] ) # radius (km)
            VeRaProfileOriginal [1][iLevel] = float ( fileLineElements [14] ) # Temperature (K)
            VeRaProfileOriginal [2][iLevel] = float ( fileLineElements [15] ) # 1-sigma Temperature uncertainty (K)
            VeRaProfileOriginal [3][iLevel] = float ( fileLineElements [8] )  * PascalToBar # pressure (bar)
            VeRaProfileOriginal [4][iLevel] = float ( fileLineElements [9] ) # 1-sigma pressure uncertainty (bar)
            VeRaProfileOriginal [5][iLevel] = float ( fileLineElements [4] ) # latitude (˚)
            VeRaProfileOriginal [6][iLevel] = float ( fileLineElements [5] ) # longitude (˚)


        fileContentTXT = HandyTools.getTextFileContent (VeRaTABFileName.split ('TAB')[0] + 'TXT')

        iLine = 0
        while 'values of 1 bar level' not in fileContentTXT [iLine]:
        
            iLine += 1
        

        iLine += 1
        

        # Return the results.
        return VeRaProfileOriginal, numberOfLevels

        
    
    # Calculate an average filtered VeRa profile at fixed altitudes.
    @staticmethod
    def getFilteredVeRaProfile (VeRaTABFileName, startAltitude = 6098., endAltitude = 6154., filteredAltitudeLevelsStep = 1., ):
        '''
        :param VeRaTABFileName: File name (and path) of the VeRa TAB file.
        :type VeRaTABFileName: str

        :param startAltitude: Start altitude (km) of the filtered profile, default = 6098km.
        :type startAltitude: float

        :param endAltitude: End altitude (km) of the filtered profile, default = 6154km.
        :type endAltitude: float

        :param filteredAltitudeLevelsStep: the step-size (km) of the filtered profile, default = 1km.
        :type filteredAltitudeLevelsStep: float
        
        :return: VeRaProfileFiltered, numberOfFilteredLevels, VeRaProfileOriginal, numberOfOriginalLevels
        :rtype: list (10 lists), int, list (7 lists), int


        **Description:**
        Calculate an average filtered VeRa profile at two altitudes between radius the  startAltitude  and  endAltitude with step  filteredAltitudeLevelsStep. 
        The default values result in 56 levels between 6089km and 6154km.

        Take a window with a total width of  filteredAltitudeLevelsStep  around each altitude level and average the temperatures and latitudes.

            | VeRaProfileFiltered [0] radius (km)
            | VeRaProfileFiltered [1] Temperature (K)
            | VeRaProfileFiltered [2] 1-sigma Temperature uncertainty (K)
            | VeRaProfileFiltered [3] pressure (bar)
            | VeRaProfileFiltered [4] 1-sigma pressure uncertainty (bar)
            | VeRaProfileFiltered [5] latitude (˚)
            | VeRaProfileFiltered [6] longitude (˚)
            | VeRaProfileFiltered [7] dT/dz (K/km)
            | VeRaProfileFiltered [8] 1-sigma dT/dz (K/km)
            | VeRaProfileFiltered [9] number of original levels used to calculate the average / filtered values

        '''

        # Read the content of the VeRa file, using the  VeRaTools.readVeRaTAB  method.
        VeRaProfileOriginal, numberOfOriginalLevels = VeRaTools.readVeRaTAB (VeRaTABFileName)

        # Make sure the user selects a valid start and end altitude, as well as altitude step.
        if filteredAltitudeLevelsStep <= 0:
        
            print ()
            print (' WARNING: filteredAltitudeLevelsStep has to be larger than 0km.')
            
            return None, None

        if (endAltitude - startAltitude) < filteredAltitudeLevelsStep:
        
            print ()
            print ( ' WARNING: difference between start and end altitudes needs to be larger than filteredAltitudeLevelsStep {}km.'.format (filteredAltitudeLevelsStep) )
            
            return None, None
            
            
        # Calculate altitude levels evenly spaced in km and average inside each altitude bin.
        numberOfFilteredLevels = int ( (endAltitude - startAltitude) / filteredAltitudeLevelsStep )

        # Create the NumPy array that will contain the resulting profiles.
        VeRaProfileFiltered = np.zeros ( (10, numberOfFilteredLevels) )
                
        # Calculate the values for each level.      
        for iFilteredAltitudeLevel in range (numberOfFilteredLevels):
                       
            VeRaProfileFiltered [0][iFilteredAltitudeLevel] = startAltitude + iFilteredAltitudeLevel
            
            iLevelsToAverage = \
             np.where ( np.logical_and ( VeRaProfileOriginal [0] > VeRaProfileFiltered [0][iFilteredAltitudeLevel] - ( filteredAltitudeLevelsStep / 2 ),
                                         VeRaProfileOriginal [0] < VeRaProfileFiltered [0][iFilteredAltitudeLevel] + ( filteredAltitudeLevelsStep / 2 ) ) )[0]


            # The number of levels in the filtered profiles.
            VeRaProfileFiltered [9][iFilteredAltitudeLevel] = len (iLevelsToAverage)

            # Loop over the filtered levels.
            if VeRaProfileFiltered [9][iFilteredAltitudeLevel] > 1:              

                # Average temperature and pressure and their uncertainty estimates (standard deviation).
                VeRaProfileFiltered [1][iFilteredAltitudeLevel], VeRaProfileFiltered [2][iFilteredAltitudeLevel], variance = \
                    DataTools.getAverageVarAndSDPYtoCPP ( VeRaProfileOriginal [1][iLevelsToAverage] )            
                VeRaProfileFiltered [3][iFilteredAltitudeLevel], VeRaProfileFiltered [4][iFilteredAltitudeLevel], variance = \
                    DataTools.getAverageVarAndSDPYtoCPP ( VeRaProfileOriginal [3][iLevelsToAverage] )
        
        
                # Average latitude and longitude.
                VeRaProfileFiltered [5][iFilteredAltitudeLevel] = DataTools.getAverageVarAndSDPYtoCPP ( VeRaProfileOriginal [5][iLevelsToAverage] )[0]
                VeRaProfileFiltered [6][iFilteredAltitudeLevel] = DataTools.getAverageVarAndSDPYtoCPP ( VeRaProfileOriginal [6][iLevelsToAverage] )[0]

            
            # In case there is only one level within the altitude bin, take the value corresponding to this level.
            elif VeRaProfileFiltered [9][iFilteredAltitudeLevel] == 1:
            
                VeRaProfileFiltered [1][iFilteredAltitudeLevel] = VeRaProfileOriginal [1][iLevelsToAverage]
                VeRaProfileFiltered [2][iFilteredAltitudeLevel] = np.nan

                VeRaProfileFiltered [3][iFilteredAltitudeLevel] = VeRaProfileOriginal [3][iLevelsToAverage]
                VeRaProfileFiltered [4][iFilteredAltitudeLevel] = np.nan
 
                VeRaProfileFiltered [5][iFilteredAltitudeLevel] = VeRaProfileOriginal [5][iLevelsToAverage]
                VeRaProfileFiltered [6][iFilteredAltitudeLevel] = VeRaProfileOriginal [6][iLevelsToAverage]
            
            # In case there are no levels in the altitude bin, set all values to NaN.
            else:

                VeRaProfileFiltered [1][iFilteredAltitudeLevel] = np.nan
                VeRaProfileFiltered [2][iFilteredAltitudeLevel] = np.nan

                VeRaProfileFiltered [3][iFilteredAltitudeLevel] = np.nan
                VeRaProfileFiltered [4][iFilteredAltitudeLevel] = np.nan
 
                VeRaProfileFiltered [5][iFilteredAltitudeLevel] = np.nan
                VeRaProfileFiltered [6][iFilteredAltitudeLevel] = np.nan
            
                                
            
        # dT/dz [7] and corresponding uncertainty estimate [8].           
        for iFilteredAltitudeLevel in range (numberOfFilteredLevels - 1):
        
            VeRaProfileFiltered [7][iFilteredAltitudeLevel] = \
             ( VeRaProfileFiltered [1][iFilteredAltitudeLevel + 1] - VeRaProfileFiltered [1][iFilteredAltitudeLevel] ) / filteredAltitudeLevelsStep

            VeRaProfileFiltered [8][iFilteredAltitudeLevel] = \
             np.sqrt ( VeRaProfileFiltered [2][iFilteredAltitudeLevel + 1]**2 + VeRaProfileFiltered [2][iFilteredAltitudeLevel]**2 ) / filteredAltitudeLevelsStep
             
        
        VeRaProfileFiltered [7][numberOfFilteredLevels - 1] = VeRaProfileFiltered [7][numberOfFilteredLevels - 2]


        return VeRaProfileFiltered, numberOfFilteredLevels, VeRaProfileOriginal, numberOfOriginalLevels
        


    # Create the table with the name  tableFileName  from the list of .TXT files present in the  topDirectory .
    @staticmethod
    def createVeRaProfilesTable (VeRaTableFileName, topDirectory, extension = 'TXT'):
        '''
        :param VeRaTableFileName: file name for the table to be written.
        :type VeRaTableFileName: str

        :param topDirectory: path name of the top directory in which the .TXT files are stored.
        :type topDirectory: str
        
        :param extension: extension of the .TXT files, default = 'TXT'
        :type extension: str


        **Description**: Create the tables with the orbit IDs and information at the one-bar pressure level of the corresponding VeRa 
        temperature profiles. The information is the Day Of Year, time of observation, Local Solar Time, Latitude, Longitude, and Solar 
        Zenith Angle at the one-bar pressure level.
        '''
    
        # Retrieve the list of .TXT files names in the sub-directories of the  topDirectory .
        listOfVeRaProfileTXTList = HandyTools.getFilesInDirectoryTree (topDirectory, extension = extension)
    
        # Go through the list and create the table.
        if listOfVeRaProfileTXTList:
        
    
            fileOpen = open (VeRaTableFileName, 'w')
            
            print (' ', file = fileOpen)
            print (' File: {}'.format ( VeRaTableFileName.replace ('../', '') ), file = fileOpen)
            print (' Created at {}'.format ( HandyTools.getDateAndTimeString () ), file = fileOpen)
            print (' ', file = fileOpen)
            print ('  Include {} files from: {}'.format (extension, topDirectory), file = fileOpen)
            print (' ', file = fileOpen)
            print ('  Values are for 1-bar level, surface radius of Venus is {:6.1f}km'.format (venusSurfaceRadius), file = fileOpen)
            print (' ', file = fileOpen)
            print (' Orbit          DOY               UTC            Altitude   Temperature   Local Solar Time  Latitude   Longitude    Solar Zenith Angle', file = fileOpen)
            print ('             yyyy-mm-dd        hh:mm:ss            (km)         (K)             (h)           (˚)         (˚)           (˚)', file = fileOpen)
            print ('C_END', file = fileOpen)
            
                
            for listOfVeRaProfileTXT in sorted (listOfVeRaProfileTXTList):
            
                content = HandyTools.getTextFileContent (listOfVeRaProfileTXT)
                
                # The .TXT files have a section starting with the text  'values of 1 bar level'  that contain the information extracted here.
                iLine = 0
                while 'values of 1 bar level' not in content [iLine]:
    
                    iLine += 1
    
                iLine += 2
                
                measurementTime = content [iLine].split (':  ')[-1]
                
                iLine += 2
                latitudeOneBar = float ( content [iLine].split (':')[-1] )
    
                iLine += 1
                longitudeOneBar = float ( content [iLine].split (':')[-1] )
 
                iLine += 1
                radiusOneBar = float ( content [iLine].split (':')[-1] )
                
                iLine += 2
                temperatureOneBar = float ( content [iLine].split (':')[-1] )
                    
                iLine += 2
                localTrueSolarTime = float ( content [iLine].split (':')[-1] )
                
                iLine += 1
                solarZenithAngle = float ( content [iLine].split (':')[-1] )
                
                
                print ( '  {}       {}        {}        {:4.1f}       {:6.2f}           {:5.2f}        {:6.2f}      {:6.2f}        {:6.2f}'.
                 format ( listOfVeRaProfileTXT.split ('/') [-2][3:7],
                          measurementTime.split ('T')[0],
                          measurementTime.split ('T')[-1],
                          radiusOneBar - venusSurfaceRadius,
                          temperatureOneBar,
                          localTrueSolarTime, 
                          latitudeOneBar,
                          longitudeOneBar,
                          solarZenithAngle ), file = fileOpen )                
                    
            
            fileOpen.close ()
           
           
        else:
        
            print ()
            print ('WARNING: no valid files found')



    # Read the values from the VeRa tables as created by the  createVeRaProfilesTable  method.
    @staticmethod
    def readValuesFromVeRaTable (VeRaTableFileName):
        '''
        :param VeRaTableFileName: file name for the table to be read.
        :type VeRaTableFileName: str

        **Description**: Read a table as created by the  createVeRaProfilesTable  of this static class and return the values of each column as lists.
        Time of day is converted to hours (= hour + minutes / 60 + seconds / 3600).
        '''
    
        VeRaTableContent = HandyTools.getTextFileContent (VeRaTableFileName)
    
        iLine = 0
        while 'C_END' not in VeRaTableContent [iLine] and iLine < len (VeRaTableContent):
        
            iLine += 1
            
        iLine += 1
        

        # Initialise empty arrays.  
        orbitID = []
        dayOfYear = []
        timeOfDay = []
        localSolarTime = [] 
        altitude = []
        temperature = []   
        longitude = []
        latitude = []
        solarZenithAngle = []
        while iLine < len (VeRaTableContent):
                
            elements = [ stringElement  for stringElement in VeRaTableContent [iLine].split (' ')  if stringElement != '' ]
    
            orbitID.append ( elements [0] ) 
            dayOfYear.append ( elements [1] )
            timeOfDay.append ( float ( elements [2].split (':')[0] ) + 
                               float ( elements [2].split (':')[1] ) / 60 + 
                               float ( elements [2].split (':')[2] ) / 3600 )
            altitude.append ( float ( elements [3] ) )
            temperature.append ( float ( elements [4] ) )
            localSolarTime.append ( float ( elements [5] ) )
            latitude.append ( float ( elements [6] ) )
            longitude.append ( float ( elements [7] ) )
            solarZenithAngle.append ( float ( elements [8] ) )

            iLine += 1
    
    
        return orbitID, dayOfYear, timeOfDay, localSolarTime, altitude, temperature, latitude, longitude, solarZenithAngle
  


        
        