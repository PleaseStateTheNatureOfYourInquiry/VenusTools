# Author: Maarten Roos-Serote
# ORCID author: 0000 0001 5001 1347

# Version: v20240604


# Standard imports.
import numpy as np
import copy

# Custom imports.
from HandyTools import HandyTools
from DataTools import DataTools

# The  planetaryimage  module can be found at https://planetaryimage.readthedocs.io/en/latest/index.html .
from planetaryimage import PDS3Image



# This is a Python class to wrangle Venus Express VMC data.
class VMCTools:
    '''
    This is a Python class to wrangle Venus Express VMC data.
    '''
    

    # Read the content of a VMC image and geocube.
    @staticmethod
    def readVMCImageAndGeoCube (VMCImageFileName):
        '''
        :param VMCImageFileName: file name (and path) of the VMC image file. It is assumed the .IMG and .GEO have the same file name.
        :type VMCImageFileName: str

        :return: VMC image as PDS3Image and flattened NumPy array, and the VMC geocube as PDS3Image as list of five flattened NumPy arrays.
        :rtype: planetaryimage.pds3image.PDS3Image, NumPy array, planetaryimage.pds3image.PDS3Image, list [NumPy array x 5]
        
        **Description:**
        Read the content of a VMC image (.IMG) and corresponding geo-cube (.GEO).
        
        The VMCImage is a PDS3Image object (`planetaryimage module <https://planetaryimage.readthedocs.io/en/latest/index.html>`_), with several attributes.
        The *image* attribute contains the uncalibrated image. For more easy manipulation, 
        the image is also returned as flattened NumPy array.
        
        The VMCGeoCube is a PDS3Image object, with several attributes. The *data* attribute contains the planes:
        
            | VMCGeoCube.data [0] Incidence Angle
            | VMCGeoCube.data [1] Emission Angle
            | VMCGeoCube.data [2] Phase Angle
            | VMCGeoCube.data [3] Latitude
            | VMCGeoCube.data [4] Longitude
        
        For more easy manipulation, each plane is returned as a flattened NumPy array, and the five arrays are packed in the  VMCGeoArraysFlattened  list.                
        '''
        
        # Make sure to not have the extension attached to the file name.
        if '.IMG' in VMCImageFileName or '.GEO' in VMCImageFileName:
        
            VMCImageFileName = VMCImageFileName [:-4]


        # Open the VMC .IMG file and create the         
        VMCImage = PDS3Image.open (VMCImageFileName + '.IMG')
        VMCImageFlattened = VMCImage.image.flatten ()
        
        
        VMCGeoCube = PDS3Image.open (VMCImageFileName + '.GEO')  
                     
        VMCGeoArraysFlattened  = []
        VMCGeoArraysFlattened.append ( VMCGeoCube.data [0].flatten () )
        VMCGeoArraysFlattened.append ( VMCGeoCube.data [1].flatten () )
        VMCGeoArraysFlattened.append ( VMCGeoCube.data [2].flatten () )      
        VMCGeoArraysFlattened.append ( VMCGeoCube.data [3].flatten () )
        VMCGeoArraysFlattened.append ( VMCGeoCube.data [4].flatten () )

        # The longitudes must run from 0˚ through 360˚ and not -180˚ through -180˚. 
        iNegativeLongitudes = np.where ( np.logical_and (VMCGeoArraysFlattened [4] >= -180, VMCGeoArraysFlattened [4] < 0) ) [0]
        VMCGeoArraysFlattened [4][iNegativeLongitudes] += 360
        VMCGeoCube.data [4] = VMCGeoArraysFlattened [4].reshape (512,512)
                
        return VMCImage, VMCImageFlattened, VMCGeoCube, VMCGeoArraysFlattened

        

    # Read the content of a VeRa .TAB file and return a numpy array.
    @staticmethod
    def VMCPhotometry ( VMCImage, 
                        VMCImageFlattened,
                        VMCGeoCube,
                        VMCGeoArraysFlattened,
                        incidenceAngleLimit = 89,
                        emissionAngleLimit = 89,
                        applyLambertLaw = True,
                        silent = False ):
        '''
        :param VMCImage: PDS3Image object (`planetaryimage module <https://planetaryimage.readthedocs.io/en/latest/index.html>`_) from reading a VMC .IMG file.
        :type VMCImage: planetaryimage.pds3image.PDS3Image

        :param VMCImageFlattened: flattened 1D array of the calibrated radiance factors of the entire VMC image.
        :type VMCImageFlattened: 1D NumPy array

        :param VMCGeoCube: PDS3Image object (`planetaryimage module <https://planetaryimage.readthedocs.io/en/latest/index.html>`_) from reading a VMC .GEO file.
        :type VMCGeoCube: planetaryimage.pds3image.PDS3Image

        :param VMCGeoArraysFlattened: flattened 1D array of the calibrated radiance factors of the entire VMC image.
        :type VMCGeoArraysFlattened: list [NumPy array x 5], see description of :py:meth:`~.readVMCImageAndGeoCube`

        :param emissionAngleLimit: valid pixels must have emission angle smaller or equal to  incidenceAngleLimit .
        :type emissionAngleLimit: float

        :param incidenceAngleLimit: valid pixels must have incidence angle smaller or equal to  incidenceAngleLimit .
        :type incidenceAngleLimit: float

        :param applyLambertLaw: apply Lambert limbdarkening law only, default = True (at this time, no other law has been programmed yet).
        :type applyLambertLaw: bool

        :param silent: print information on run, default = False.
        :type silent: bool

        :return: calibrated VMC image, flattened calibrated VMC image, average incidence angle, standard deviation incidence angle, average emission angle, standard emission incidence angle, average phase angle, standard deviation phase angle, radiance scaling factor 
        :rtype: 2D NumPy array (512x512), 1d NumPy array, float, float, float, float, float, float,
        
        **Description:**
        This method is used to calibrate a VMC image.
        The calibrated radiance factor :math:`RF_{x,y}` for a valid (= on Venus disk with incidence angle < :code:`incidenceAngleLimit` and emission 
        angle < :code:`emissionAngleLimit`) pixel :math:`(x,y)` in a VMC image is:
                
            :math:`RF_{x,y} = \\pi \\beta R_{observed - x,y} \\frac {d_{Venus}}{S_{Sun}}`
        
        
        where :math:`\\beta` is the calibration correction factor (see `Shalygina  et al. 2015 <https://doi.org/10.1016/j.pss.2014.11.012>`_, their Table 1), :math:`R_{observed - x,y}` is the value at the pixel in ADU times the radiance scaling factor read from the VMC image header (:code:`VMCImage.label ['RADIANCE_SCALING_FACTOR'].value`, when read with `planetaryimage module <https://planetaryimage.readthedocs.io/en/latest/index.html>`_ in Python) in :math:`W/m^2/\\mu m/ster/ADU`, :math:`d_{Venus}` is the distance of Venus to the Sun in AU and :math:`S_{Sun}` the solar flux in :math:`W/m^2/\\mu m` at 1AU (see `Lee et al. 2015 <http://dx.doi.org/10.1016/j.icarus.2015.02.015>`_ their Equation 2).
        
        
        The :math:`\\beta` factor is 2.34 for orbits before 2639 and 1 for later orbits.
        
        For :math:`S_{Sun}` (from `Lee et al. 2015 <http://dx.doi.org/10.1016/j.icarus.2015.02.015>`_ their Equation 1):
        
            :math:`S_{Sun} = \\frac {\\int S_{irradiance}(\\lambda) T(\\lambda) d\\lambda}{\\int T (\\lambda) d\\lambda}`
        
        The :math:`S_{irradiance}` is determined using the `Solar Spectra website <https://www.nrel.gov/grid/solar-resource/spectra.html>`_ and the transmission function :math:`T (\\lamda)` is taken from a parametrisation of the UV part of the transmission function of the VMC camera as published by `Markiewicz et al. 2007 <https://doi.org/10.1016/j.pss.2007.01.004>`_ their Figure 3.
        
        This results in a value of :math:`S_{Sun} = 1081 W/m^2/\\mu m`.
        
        The value of :math:`d_{Venus} = 0.723AU` to within 1% in :math:`d_{Venus}^2` over the orbit of Venus.
        
        Finally, a Lambert limb-darkening law is applied. 
        
        The calibration is only done for pixels that fall on the Venus disk, which have valid longitudes attached to them and an incidence angle of less than or equal to the limit set by the user (default = 89˚).
        
        The first four input variables come from reading a VMC image with the :py:meth:`~.readVMCImageAndGeoCube` method of this class.
        '''
    
        # Make sure the longitudes are running from 0 through 360˚ and not -180˚ through 180˚.
        iLongitude = np.where ( np.logical_and ( VMCGeoArraysFlattened [4] >= -180, VMCGeoArraysFlattened [4] < 0 ) ) [0]
        VMCGeoArraysFlattened [4][iLongitude] += 360
        
        # Collect the indices of all the points on the illuminated part of the Venus disk.
        iOnDiskValid1 = np.where ( np.logical_and ( abs ( VMCGeoArraysFlattened [4] ) <= 360, 
                                                    abs ( VMCGeoArraysFlattened [0] ) <= incidenceAngleLimit ) ) [0]

        
        iOnDiskValid2 = np.where ( VMCGeoArraysFlattened [1][iOnDiskValid1] <= emissionAngleLimit ) [0]
                                                    
        # Valid points with incidence and emission angle limits.
        iOnDiskValid = iOnDiskValid1 [iOnDiskValid2]
        
   
        # Extract the Radiance Scaling Factor in the right units of W/m2/ster/micron, instead of the reported units W/m3/ster.   
        radianceScalingFactor = VMCImage.label ['RADIANCE_SCALING_FACTOR'].value / 1000000.        

        # The beta factor is ...
        betaFactor = 2.34  if VMCImage.label ['ORBIT_NUMBER'] <= 2638  else  1

        VMCImageCalibratedFlattened = np.zeros ( len (VMCImageFlattened) ) - 1.
        VMCImageCalibratedFlattened [iOnDiskValid] = radianceScalingFactor * VMCImageFlattened [iOnDiskValid] * betaFactor * np.pi * (0.723 * 0.723) / 1081


        # Use Lambert's law for the correction of the indicence angle.
        if applyLambertLaw:
     
            VMCImageCalibratedFlattened [iOnDiskValid] /= np.cos  ( VMCGeoArraysFlattened [0][iOnDiskValid] * np.pi / 180 )


        if not silent:
        
            print ()
            print ( ' Calibrating image {}'.format (VMCImage.filename) )
            print ( '  - valid point when longitude 0˚ - 360, incidence angle < {}˚ and emission angle < {}˚;'.format (incidenceAngleLimit, emissionAngleLimit) )
            print ( '  - number of valid points = {:6d};'.format ( len (iOnDiskValid) ) )
            print ( '  - radiance scaling factor = {} W/m2/micron/ster;'.format (radianceScalingFactor) )

            if applyLambertLaw:

                print ( '  - applying Lambert\'s law.')


        incidenceAngeAverage = DataTools.getAverageVarAndSDPYtoCPP ( VMCGeoArraysFlattened [0][iOnDiskValid] )
        emissionAngleAverage = DataTools.getAverageVarAndSDPYtoCPP ( VMCGeoArraysFlattened [1][iOnDiskValid] )
        phaseAngleAverage = DataTools.getAverageVarAndSDPYtoCPP ( VMCGeoArraysFlattened [2][iOnDiskValid] )

        return VMCImageCalibratedFlattened.reshape (512, 512), VMCImageCalibratedFlattened, \
               incidenceAngeAverage [0], incidenceAngeAverage [1], \
               emissionAngleAverage [0], emissionAngleAverage [1], \
               phaseAngleAverage [0], phaseAngleAverage [1], radianceScalingFactor



    # 
    @staticmethod
    def getWindAdvectedBox (latitudeVeRaSounding, longitudeVeRaSounding, timeDifferenceHours, oneSigmaZonalWind = 20, oneSigmaMeridionalWind = 12):
        '''
        :param latitudeVeRaSounding: latitude (˚) of the VeRa sounding location.
        :type latitudeVeRaSounding: float

        :param longitudeVeRaSounding: longitude (˚) of the VeRa sounding location.
        :type longitudeVeRaSounding: float

        :param timeDifferenceHours: difference in hours between the VeRa sounding and the VMC image acquisition.
        :type timeDifferenceHours: float        

        :param oneSigmaZonalWind: the standard deviation for the zonal wind. Default is 20m/s from Khatunstev et al. (2013) Figure 10a.
        :type oneSigmaZonalWind: float        

        :param oneSigmaMeridionalWind: the standard deviation for the meridional wind. Default is 12m/s from Khatunstev et al. (2013) Figure 10b.
        :type oneSigmaMeridionalWind: float        


        :return: latitudeCentre, latitudeLimits, longitudeCentre, longitudeLimits.
        :rtype: float, list [float, float], float, list [ [float, float], [float, float] ]
        
        **Description**:   
        Given the latitude and longitude of the VeRa sounding, as well as the time difference between the VeRa sounding and the VMC image, determine the
        latitude-longitude box where the atmosphere of the VeRa sounding was positioned at the time of the VMC image recording.
        
        The parametrisation is based on `Khatuntsev et al. (2013) <http://dx.doi.org/10.1016/j.icarus.2013.05.018s>`_ figures 10a (zonal) and 10b (meridional).
        The uncertainties (standard deviation) `oneSigmaZonalWind` and `oneSigmaMeridionalWind` are estimated based on the grey areas in those figures.
        '''
                        

        # conversion factor for metres to degrees longitude at the equator at 70km altitude (radius of Venus taken to be 6052km). 
        degreesLatitudePerMetre = 180 / ( np.pi * (6052 + 70) * 1000 )
        degreesLongitudePerMetreAtVeRaLatitude = degreesLatitudePerMetre / np.cos ( np.pi * latitudeVeRaSounding / 180 )
        timeDifferenceSeconds = timeDifferenceHours * 3600
            
        if latitudeVeRaSounding <= -50:
        
            longitudeCentre = longitudeVeRaSounding + degreesLongitudePerMetreAtVeRaLatitude * timeDifferenceSeconds * ( -94 - (latitudeVeRaSounding + 50) * (65.6 / 25) )
                
        elif latitudeVeRaSounding <= -40 and latitudeVeRaSounding >= -50:
        
            longitudeCentre = longitudeVeRaSounding + degreesLongitudePerMetreAtVeRaLatitude * timeDifferenceSeconds * ( -101.5 - (latitudeVeRaSounding + 40) * (7.5 / 10) )
        
        elif latitudeVeRaSounding <= -15 and latitudeVeRaSounding >= -40:
        
            longitudeCentre = longitudeVeRaSounding + degreesLongitudePerMetreAtVeRaLatitude * timeDifferenceSeconds * ( -93 + (latitudeVeRaSounding + 15) * (8.5 / 35) )
               
        elif latitudeVeRaSounding > -15:
        
            longitudeCentre = longitudeVeRaSounding + degreesLongitudePerMetreAtVeRaLatitude * timeDifferenceSeconds * -93
        

        longitudeHalfRange = abs (degreesLongitudePerMetreAtVeRaLatitude * timeDifferenceSeconds * oneSigmaZonalWind)
        if longitudeVeRaSounding > longitudeCentre:
        
            longitudeBorderMaximum = longitudeCentre + longitudeHalfRange
            longitudeBorderMinimum = longitudeCentre - longitudeHalfRange
            
        else:
        
            longitudeBorderMaximum = longitudeCentre + longitudeHalfRange
            longitudeBorderMinimum = longitudeCentre - longitudeHalfRange



        longitudeLimits = [ [longitudeBorderMinimum, longitudeBorderMaximum], [longitudeBorderMinimum, longitudeBorderMaximum] ]
        # All values are negative, the box lies entirely beyond the 0˚ meridian.
        if longitudeBorderMaximum < 0:

            longitudeBorderMaximum += 360
            longitudeCentre += 360
            longitudeBorderMinimum += 360
                    
            longitudeLimits = [ [longitudeBorderMinimum, longitudeBorderMaximum], [longitudeBorderMinimum, longitudeBorderMaximum] ]


        # The centre and the minimum border are beyond the 0˚ meridian.            
        if longitudeCentre < 0 and longitudeBorderMaximum > 0:
        
            longitudeCentre += 360
            longitudeBorderMinimum += 360                             

            longitudeLimits = [ [longitudeBorderMinimum, 360], [0, longitudeBorderMaximum] ]


        # Only the minimum border is beyond the 0˚ meridian. 
        if longitudeBorderMinimum < 0 and longitudeCentre > 0:
        
            longitudeBorderMinimum += 360                             

            longitudeLimits = [ [longitudeBorderMinimum, 360], [0, longitudeBorderMaximum] ]



        # All values are > 360˚, the box lies entirely beyond the 360˚ meridian.
        if longitudeBorderMinimum > 360:

            longitudeBorderMaximum -= 360
            longitudeCentre -= 360
            longitudeBorderMinimum -= 360

            longitudeLimits = [ [longitudeBorderMinimum, longitudeBorderMaximum], [longitudeBorderMinimum, longitudeBorderMaximum] ]
        

        # The centre and the maximum border are beyond the 360˚ meridian.            
        if longitudeCentre > 360 and longitudeBorderMinimum < 360:
        
            longitudeCentre -= 360
            longitudeBorderMaximum -= 360                             

            longitudeLimits = [ [longitudeBorderMinimum, 360], [0, longitudeBorderMaximum] ]


        # Only the maximum border is beyond the 360˚ meridian.
        if longitudeCentre < 360 and longitudeBorderMaximum > 360:
        
            longitudeBorderMaximum -= 360                             

            longitudeLimits = [ [longitudeBorderMinimum, 360], [0, longitudeBorderMaximum] ]



        if latitudeVeRaSounding <= -75:
        
            latitudeCentre = latitudeVeRaSounding
                
        elif latitudeVeRaSounding <= -50 and latitudeVeRaSounding > -75:
        
            latitudeCentre = latitudeVeRaSounding + degreesLatitudePerMetre * timeDifferenceSeconds * ( -9.58 - (latitudeVeRaSounding + 50) * (9.38 / 25) )
               
        elif latitudeVeRaSounding <= -20 and latitudeVeRaSounding > -50:
        
            latitudeCentre = latitudeVeRaSounding + degreesLatitudePerMetre * timeDifferenceSeconds * ( -6.5 + (latitudeVeRaSounding + 20) * (3.08 / 30) )
               
        elif latitudeVeRaSounding > -20:
        
            latitudeCentre = latitudeVeRaSounding + degreesLatitudePerMetre * timeDifferenceSeconds * ( -3.26 + latitudeVeRaSounding * (3.24 / 20) )


        latitudeHalfRange = abs (degreesLatitudePerMetre * timeDifferenceSeconds * oneSigmaMeridionalWind)

        latitudeBorderMinimum = max (-90, latitudeCentre - latitudeHalfRange)
        latitudeBorderMaximum = max (-90, latitudeCentre + latitudeHalfRange)
        
        latitudeLimits  = [latitudeBorderMinimum, latitudeBorderMaximum]

        return latitudeCentre, latitudeLimits, longitudeCentre, longitudeLimits 



    # Define the colours for each Venus Express mission section for the scatter plots.
    def getColourForVEXMissionSection (orbitOrImageName):
        '''
        :param orbitOrImageID: ID of a VMC image, for example 'V2811_0080_UV2'
        :type orbitOrImageID: str or int        

        :return: colour string for use in matplotlib.
        :rtype: str
        
        **Description:**
        Define the colours for each Venus Express mission section for making scatter plots.
        The information about the orbit IDs per mission section comes from the file 
        :file:`/Users/maarten/Science/Venus/Temperature-UV_Analysis_2024/Data/VEX/VEX-SCIOPS-LI-053_1_1_VEX_Orbit_Date_DOY_Listing_2014Sep15.numbers`,
        see also the link to `Venus Express at the ESA PSA <https://www.cosmos.esa.int/web/psa/venus-express>`_
        '''


        # Make sure  orbitID  is a string variable.
        if type (orbitOrImageName) != str:
        
            orbitOrImageName = str (orbitOrImageName)



        # Extract the VEX orbitID for the string or number provided by the user.
        if '.' in orbitOrImageName:
        
            orbitOrImageName = orbitOrImageName.split ('.')[0]


        if '_' in orbitOrImageName:
        
            orbitOrImageName = orbitOrImageName.split ('_')[0]
        
        
        orbitID = orbitOrImageName.replace ('V', '')
        

        # Determine the colour for the mission section of the image / orbit.

        # Nominal mission.
        if orbitID <= '0547':
        
           colourString = 'green'     
            
        # Extension 1.
        elif '0548' <= orbitID <= '1135':
        
           colourString = 'blue'
    
        # Extension 2.
        elif '1136' <= orbitID <= '1583':
        
           colourString = 'purple'
    
        # Extension 3.
        elif '1584' <= orbitID <= '2451':
        
           colourString = 'black'
    
        # Extension 4, including the South Polar Dynamics Campaign.
        elif '2452' <= orbitID <= '3537':
        
           colourString = 'grey'
           
           # South Polar Dynamics Campaign.
           if '2775' <= orbitID <= '2811':
        
                colourString = 'red'
           
    
        return colourString


