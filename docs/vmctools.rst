.. _vmctools:

VMCTools
====================


| :py:meth:`~.readVMCImageAndGeoCube`
| :py:meth:`~.VMCPhotometry`
| :py:meth:`~.getWindAdvectedBox`
| :py:meth:`~.getColourForVEXMissionSection`



.. automethod:: VMCTools.VMCTools.readVMCImageAndGeoCube


    The longitude values stored in the .GEO files go from -180˚ to +180˚. In this function the longitude values are transformed so that they run between 0˚ to +360˚. 
    In this plot it is shown the transformation is done correctly.
    
    .. figure:: ../andere/images/longitudeTransformationCheck.png


.. automethod:: VMCTools.VMCTools.VMCPhotometry



.. automethod:: VMCTools.VMCTools.getWindAdvectedBox


    .. admonition:: Parametrisation of the wind profiles.
    
        Khatuntsev et al. (2013) report on the zonal and meridional wind profiles measured from VEX orbits up to 2299 (10-year period) and
        present the result in their figures 10(a) and (b): *Mean zonal (a) and meridional (b) profiles of the wind speed derived over the period of 10 venusian years by manual cloud tracking. Error bars correspond to 99.9999% 5σ-x confidence interval based on the standard deviation of the weighted mean. Standard deviations are presented by shadowed areas.*
        
        They present the formulae (their Equations (1) and (2)) with which the winds have been calculated by comparing two images:
        
        .. math::
        
            U = \frac {(\lambda_2 - \lambda_1) (R+h) cos (\theta)}{\Delta t}
            
        
        .. math::
        
            V = \frac {(\theta_2 - \theta_1) (R + h)}{\Delta t}
    
        where the indices :math:`1` and :math:`2` refer to the first and second image, :math:`\lambda` is the longitude and :math:`\theta` the latitude. Since :math:`U` is negative (see table below), it means the wind blows the clouds in the direction of smaller Venus longitudes.
    
        
        From my (physical) notebook entry on 19-03-2015: the average zonal wind is determined from figure 10a in Khatuntsev et al. (2013)
        and can be parametrised as (:math:`U` in units of m/s):
            
        .. csv-table:: **Parametrised average zonal wind**
           :header: "Latitude range (˚)", "wind parametrisation (m/s)", notes
     
            "(-75˚, -50˚]", :math:`U({\theta})` =  -94   + (:math:`{\theta}` + 50) * (65.6/-25), also use for latitudes down to -90˚.
            "(-50˚, -40˚]", :math:`U({\theta})` = -101.5 + (:math:`{\theta}` + 40) * (7.5/-10)
            "(-40˚, -15˚]",  :math:`U({\theta})` =  -93 + (:math:`{\theta}` + 15) * (-8.5/-35)
            "(-15˚, 0˚]", :math:`v{\theta})` = -93
    
    
        .. figure:: //Users/maarten/Science/Venus/Temperature-UV_Analysis_2024/Analysis/VMC/Step02/KhatuntsevWindProfiles/Khatuntsev_2013_Figure10a_ZonalWind.jpg
            :scale: 5% 
        
            Khatuntsev et al. (2013) Figure 10a. Zonal wind parametrisation.   
    
    
        From my (physical) notebook entry on 24-04-2015: the average meridional wind is determined from figure 10b in Khatuntsev et al. (2013)
        and can be parametrised as (:math:`V` in units of m/s):
    
    
        .. csv-table:: **Parametrised average meridional wind**
           :header: "Latitude range (˚)", "wind parametrisation (m/s)"
    
            "(-90˚, -75˚]",  :math:`V({\theta})`) =  0
            "(-75˚, -50˚]",  :math:`V({\theta})`) = -9.58 + ( :math:`{\theta}` + 50 ) * (9.38/-25)
            "(-50˚, -20˚]",  :math:`V({\theta})`) = -6.5  + ( :math:`{\theta}` + 20 ) * (-3.08/-30)
            "(-20˚, 0˚]",    :math:`V({\theta})`) = -3.26 + ( :math:`{\theta}` + 0 )   * (-3.24/-20)
    
        .. figure:: //Users/maarten/Science/Venus/Temperature-UV_Analysis_2024/Analysis/VMC/Step02/KhatuntsevWindProfiles/Khatuntsev_2013_Figure10b_MeridionalWind.jpg
            :scale: 5% 
            
            Khatuntsev et al. (2013) Figure 10b. Meridional wind parametrisation.   
        
        


.. automethod:: VMCTools.VMCTools.getColourForVEXMissionSection




