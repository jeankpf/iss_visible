from math import radians, degrees
import requests

import orekit
VM = orekit.initVM()
from orekit.pyhelpers import setup_orekit_curdir
setup_orekit_curdir()
from org.orekit.frames import FramesFactory, TopocentricFrame
from org.orekit.bodies import OneAxisEllipsoid, GeodeticPoint, CelestialBodyFactory
from org.orekit.propagation.analytical.tle import TLE, TLEPropagator
from org.orekit.utils import IERSConventions, Constants, PVCoordinatesProvider
from org.orekit.propagation.events import ElevationDetector, EventsLogger
from org.orekit.propagation.events.handlers import ContinueOnEvent


def get_station(longi, lat, name, planet=OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                                                          Constants.WGS84_EARTH_FLATTENING,
                                                          FramesFactory.getITRF(IERSConventions.IERS_2010, True))):
    """
    Returns the wanted Topocentric Frame computed thanks to its coordinates 

    Parameters
    ----------
    longi : float
        longitude of the wanted Topocentric Frame.
    lat : float
        Latitude of the wanted Topocentric Frame.
    name : str
        Wanted name for the Topocentric Frame.
    planet : OnesAxisEllipsoid, optional
        Planet where the Topocentric Frame should be associated to. The default is our planet, the Earth.

    Returns
    -------
    station_frame : Topocentric Frame
        Wanted Topocentric Frame.

    """
    longitude = radians(longi)
    latitude = radians(lat)
    station = GeodeticPoint(latitude, longitude, 0.0)
    station_frame = TopocentricFrame(planet, station, name)
    return station_frame

# Looks for ISS TLE on internet (hopefully this page is updated quite often!)
def get_iss_tle():
    """
    Returns the Two Line Elements (TLE) of the ISS.

    Returns
    -------
    tle : TLE
        Two Line Elements of the ISS.

    """
    url = 'https://www.celestrak.com/NORAD/elements/stations.txt'
    r = requests.get(url, allow_redirects=True)
    open('stations.txt', 'wb').write(r.content)
    file = open('stations.txt', 'r')
    lines = file.readlines()
    # assert lines[0] == "ISS (ZARYA)"
    # tle1 = "1 25544U 98067A   20235.39741112  .00007958  00000-0  15147-3 0  9990"
    # tle2 = "2 25544  51.6464  20.8011 0001866  40.7927 349.0439 15.49187314242283"
    tle1 = lines[1]
    tle2 = lines[2]
    tle = TLE(tle1, tle2)
    assert isinstance(tle, TLE)
    file.close()
    return tle




def is_iss_visible(e_deg):
    """
    Arbitraty criteria on the elevation value to check if the ISS is visble.

    Parameters
    ----------
    e_deg : float
        Elevation angle in degrees.

    Returns
    -------
    bool
        True if the ISS is visible, False otherwise.

    """
    if e_deg>5:
        return True
    return False

def when_is_iss_visible(longi, lat, station_name):
    """
    Returns some parameters of the ISS at its maximum elevation.

    Parameters
    ----------
    longi : float
        Longitude of the place where you stand on Earth.
    lat : float
        Latitude of the place where you stand on Earth.
    station_name : str
        Name of this place.

    Returns
    -------
    extrap_date : AbsoluteDate
        Date and time of the event.
    elevation : float
        Maximum elevation angle in radians.
    current_x : float
        X coordinate in the Topocentric Frame.
    current_y : TYPE
        Y coordinate in the Topocentric Frame.

    """
    iss_tle = get_iss_tle()
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                             Constants.WGS84_EARTH_FLATTENING, itrf)
    station_frame = get_station(longi, lat, station_name, planet=earth)
    propagator = TLEPropagator.selectExtrapolator(iss_tle)

    current_x = 100001
    current_y = 100001
    current_ele = 0
    extrap_date = iss_tle.getDate()
    inertial_frame = FramesFactory.getEME2000()
    while not is_iss_visible(degrees(current_ele)):
        pv = propagator.getPVCoordinates(extrap_date, inertial_frame)
        pos_tmp = pv.getPosition()
        inertial2station_frame = inertial_frame.getTransformTo(station_frame, extrap_date)
        current_pos_station_frame = inertial2station_frame.transformPosition(pos_tmp)
        current_x = current_pos_station_frame.getX()
        current_y = current_pos_station_frame.getY()
        current_ele = station_frame.getElevation(current_pos_station_frame,
                                                 station_frame, extrap_date)
        extrap_date = extrap_date.shiftedBy(10.0)
    elevation = station_frame.getElevation(pos_tmp, inertial_frame,
                                           extrap_date)
    return extrap_date, elevation, current_x, current_y

def when_is_iss_visible_local_time(longi, lat, station_name, time_zone, summer_time):
    """
    Returns the date, time and elevation of the ISS when it is visible (e.g. at its maximum elevation).

    Parameters
    ----------
    longi : float
        Longitude of the place where you stand on Earth.
    lat : float
        Latitude of the place where you stand on Earth.
    station_name : str
        Name of this place.
    time_zone : int
        Time zone of the place where you stand on Earth (e.g. 1 for Paris, 0 for London...).
    summer_time : bool
        True if it's summer time, False otherwise.

    Returns
    -------
    date : AbsoluteDate
        Date and time of the event.
    elevation : float
        Maximum elevation angle in degrees.

    """
    if summer_time:  # summer time
        time_zone += 1
    date, elevation, _, _ = when_is_iss_visible(longi, lat,
                                              station_name)
    date = date.shiftedBy(time_zone * 60.0 * 60)     # get local time
    return date, degrees(elevation)


def get_iss_ground_plot(start_date, end_date):
    return 0

def is_day_utc(utc_date, station_frame):
    """
    Check if it's day or not at a given date (in UTC time) at a given place on Earth'

    Parameters
    ----------
    utc_date : AbsoluteDate
        date and time in the UTC.
    station_frame : Topocentric Frame
        Frame of the place on Earth.

    Returns
    -------
    bool
        True if it is day, False otherwise.

    """
    j2000 = FramesFactory.getEME2000()
    pv_sun = CelestialBodyFactory.getSun()
    pv_sun = PVCoordinatesProvider.cast_(pv_sun)
    sun_pos_j2000 = pv_sun.getPVCoordinates(utc_date, j2000).getPosition()
    j2000_to_topocentric = j2000.getTransformTo(station_frame, utc_date)
    sun_pos_topo = j2000_to_topocentric.transformPosition(sun_pos_j2000)
    return sun_pos_topo.getZ() > 0

def is_day_local_time(local_date, station_frame, time_zone, summer_time):
    """
    Check if it's day or not at a given date (in local time) at a given place on Earth'


    Parameters
    ----------
    utc_date : AbsoluteDate
        date and time in the UTC.
    station_frame : Topocentric Frame
        Frame of the place on Earth.
    time_zone : int
        Time zone of the place where you stand on Earth (e.g. 1 for Paris, 0 for London...).
    summer_time : bool
        True if it's summer time, False otherwise.

    Returns
    -------
    bool
        True if it is day, False otherwise.

    """
    j2000 = FramesFactory.getEME2000()
    pv_sun = CelestialBodyFactory.getSun()
    pv_sun = PVCoordinatesProvider.cast_(pv_sun)
    if summer_time:  # heure d'ete
        time_zone += 1
    sun_pos_j2000 = pv_sun.getPVCoordinates(local_date.shiftedBy(-time_zone * 60.0 * 60),
                                            j2000).getPosition()
    j2000_to_topocentric = j2000.getTransformTo(station_frame,
                                                local_date.shiftedBy(-time_zone * 60.0 * 60))
    sun_pos_topo = j2000_to_topocentric.transformPosition(sun_pos_j2000)
    return sun_pos_topo.getZ() > 0

def when_is_iss_visible_local_time_bis(longi, lat, station_name, time_zone, summer_time):
    if summer_time:  # heure d'ete
        time_zone += 1
    # get UCT
    date, elevation = when_is_iss_visible_bis(longi, lat, station_name)
    # get local time
    date = date.shiftedBy(time_zone * 60.0 * 60)
    return date, degrees(elevation)

def when_is_iss_visible_bis(longi, lat, station_name):
    iss_tle = get_iss_tle()
    itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, True)
    earth = OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                             Constants.WGS84_EARTH_FLATTENING, itrf)
    station_frame = get_station(longi, lat, earth, station_name)
    propagator = TLEPropagator.selectExtrapolator(iss_tle)
    ele_detector = ElevationDetector(60.0, 0.001, station_frame).withConstantElevation(radians(5.0)).withHandler(ContinueOnEvent())
    logger = EventsLogger()
    logged_detector = logger.monitorDetector(ele_detector)
    # propagator = Propagator.cast_(propagator)
    propagator.addEventDetector(logged_detector)
    initial_date = iss_tle.getDate()
    propagator.propagate(initial_date, initial_date.shiftedBy(3600.0*0.5)) ##20 days here
    date_ele_max = logger.getLoggedEvents().get(0).getState().getDate()
    pos_max_j2000 = logger.getLoggedEvents().get(0).getState().getPVCoordinates().getPosition()
    ele_max = station_frame.getElevation(pos_max_j2000,
                                         FramesFactory.getEME2000(), date_ele_max)
    return date_ele_max, ele_max
