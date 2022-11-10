"""
Astrological calculations for primary directions
"""
from datetime import datetime
from math import tan, sin, cos, pi, atan2, asin
from typing import Optional
from components.vector import coords as crd
from sphere import Sphere, true_distance


def possible_aspects(aspect: Optional[int] = None):
    """
    Returns the list of possible aspects (sinister and
    dexter aspects) based on aspect's angle in degrees
    """
    if aspect is None:
        return [-120, -90, - 60, 0, 60, 90, 120, 180]
    if aspect in [60, 90, 120]:
        return [-1 * aspect, aspect]
    if aspect == 180:
        return [180]
    return [0]


class Directions:
    """
    Calculates different types of primary directions
    """

    def __init__(self, sphere: Sphere) -> None:
        self.sphere = sphere

    def quadrant(self, right_asc: float, dec: float) -> int:
        """
        Returns the celestial sphere quadrant,
        which a given point belongs to.
        0 - east-top, or east horizon, or MC
        1 - east-bottom
        2 - west-bottom
        3 - west-top, or west horizon
        """
        xyz_hrz = self.sphere.set_equatorial(right_asc, dec).horizontal_xyz()
        if xyz_hrz.x >= 0 and xyz_hrz.z >= 0:
            return 0
        if xyz_hrz.x > 0 and xyz_hrz.z < 0:
            return 1
        if xyz_hrz.x < 0 and xyz_hrz.z < 0:
            return 2
        return 3

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Placidus direction common methods:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def md_portion(self, right_asc: float, dec: float) -> Optional[float]:
        """
        Returns the meridian distance portion,
        MDP of upper/lower meridian distance
        relative to diurnal/nocturnal semiarc
        for a given point
        """
        dsa = self.sphere.dsa(dec)
        umd = self.sphere.umd(right_asc)
        if dsa is None:
            return None

        if umd <= dsa:
            return umd / dsa
        return (180 - umd) / (180 - dsa)

    def __ra_conj_placidus(self,
                           mdp: float,
                           quadrant: int,
                           dec: float) -> float:
        """
        Returns the RA of point A with
        declination dec(A), which `conjuncts`
        point B with the given meridian
        distance portion MDP(B) in a given
        quadrant(B) in placidus system
        """
        ramc = self.sphere.ramc
        dsa = self.sphere.dsa(dec)
        if dsa is None:
            return None
        nsa = 180 - dsa

        if quadrant == 0:
            return ramc + mdp * dsa
        if quadrant == 3:
            return ramc - mdp * dsa
        if quadrant == 1:
            return ramc + 180 - mdp * nsa
        return ramc + 180 + mdp * nsa

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Placidus mundane directions:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def aspect_positions_placidus_mundane(self,
                                          acceptor: crd.Vector,
                                          aspect: Optional[int] = None) -> Optional[list[float]]:
        """
        Returns mundane aspect positions of the acceptor.
        It is the point on equatorial plane with the same
        meridian distance portion in the same quadrant
        as the acceptor, plus all major aspects to
        that point on the equatorial plane.
        """

        acceptor_rasc = acceptor.equatorial().rasc
        acceptor_dec = acceptor.equatorial().dec

        acceptor_quadrant = self.quadrant(acceptor_rasc, acceptor_dec)
        acceptor_mdp = self.md_portion(acceptor_rasc, acceptor_dec)
        if acceptor_mdp is None:
            return None

        # Find RA of the equator plane with the same
        # meridian distance portion as the acceptor
        eqt_rasc = self.__ra_conj_placidus(
            acceptor_mdp, acceptor_quadrant, dec=0
        )

        # Find aspect positions of that point on the equator plane
        return [
            dict(rasc=(eqt_rasc + _a) % 360, aspect=abs(_a))
            for _a in possible_aspects(aspect)
        ]

    def placidus_mundane(self,
                         promissor: crd.Vector,
                         acceptor: crd.Vector,
                         aspect: Optional[int] = None) -> Optional[list[dict]]:
        """
        Returns mundane primary directions of the promissor
        to aspect points of the acceptor
        """
        aspect_positions = self.aspect_positions_placidus_mundane(
            acceptor, aspect)
        if aspect_positions is None:
            return None

        # Find conjunctions of the promissor with these aspect positions
        promissor_dec = promissor.equatorial().dec
        promissor_rasc = promissor.equatorial().rasc
        directions = []
        for item in aspect_positions:
            eqt_quadrant = self.quadrant(item['rasc'], dec=0)
            eqt_mdp = self.md_portion(item['rasc'], dec=0)
            aspect_rasc = self.__ra_conj_placidus(
                eqt_mdp, eqt_quadrant, promissor_dec
            )

            # In case of no ascension
            if aspect_rasc is None:
                continue
            distance = true_distance(aspect_rasc, promissor_rasc)
            check = true_distance(aspect_rasc, promissor_rasc - distance)

            # Negative distances for converse direction
            directions.append(dict(
                dist=distance if abs(check) < 1e-10 else -distance,
                aspect=item['aspect']
            ))

        return directions

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Placidus zodiacal directions:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    def aspect_positions_placidus_zodiac(self,
                                         promissor: crd.Vector,
                                         aspect: Optional[int] = None) -> Optional[list[dict]]:
        """
        Returns zodiacal aspect positions of the
        promissor. It is the projection of the
        promissor in equatorial plane, plus all
        major aspects to that point on the equatorial plane.
        """
        prom_lon = promissor.ecliptical().lon
        return [dict(lon=(prom_lon + _a) % 360, aspect=abs(_a)) for _a in possible_aspects(aspect)]

    def placidus_zodiac(self,
                        promissor: crd.Vector,
                        acceptor: crd.Vector,
                        aspect: Optional[int] = None,
                        field_plane_lat: Optional[float] = 0,
                        ) -> list[dict]:
        """
        Returns zodiacal primary directions of the promissor's
        aspect points to the acceptor
        """
        aspect_positions = self.aspect_positions_placidus_zodiac(
            promissor, aspect)

        # Find conjunctions of the promissor with these aspect positions
        directions = []
        for item in aspect_positions:
            new_promissor = self.sphere.set_ecliptical(
                item['lon'], field_plane_lat)
            arc = self.placidus_mundane(
                new_promissor, acceptor, aspect=0
            )

            # In case of no ascension:
            if not arc:
                continue
            directions.append(dict(dist=arc[0]['dist'], aspect=item['aspect']))

        return directions

    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
    # Regiomontanus mundane directions:
    # ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌

    def __oblique_asc_utp__(self, rasc: float, dec: float) -> float:
        """
        Returns oblique ascension under the pole
        for a given point.
        """
        tan_p = self.sphere.__constants__.tan_p
        tan_d = tan(dec * pi / 180)

        oa_asc = self.sphere.ramc + 90
        cos_oa = cos(oa_asc * pi / 180)
        sin_oa = sin(oa_asc * pi / 180)

        _y = sin(rasc * pi / 180) - tan_d * tan_p * cos_oa
        _x = cos(rasc * pi / 180) + tan_d * tan_p * sin_oa
        return atan2(_y, _x) * 180 / pi

    def __ra_from_oa_utp__(self, dec: float, oa_utp: float) -> float:
        """
        Returns right ascension of the point with a
        given declination and oblique ascension under the pole
        """
        tan_d = tan(dec * pi / 180)
        tan_p = self.sphere.__constants__.tan_p
        sin_t = sin((oa_utp - self.sphere.ramc) * pi / 180)

        return asin(
            tan_d * tan_p * sin_t
        ) * 180 / pi + oa_utp

    def aspect_positions_regio_mundane(self,
                                       acceptor: crd.Vector,
                                       aspect: Optional[int] = None) -> Optional[list[float]]:
        """
        Returns mundane aspect positions of the acceptor.
        It is the oblique ascension under the pole on
        equator, plus all major aspects to that point on
        the equatorial plane.
        """
        oa_utp = self.__oblique_asc_utp__(
            acceptor.equatorial().rasc,
            acceptor.equatorial().dec
        )
        # Find aspect positions of that point on the equator plane
        return [
            dict(rasc=(oa_utp + _a) % 360, aspect=abs(_a))
            for _a in possible_aspects(aspect)
        ]

    def regio_mundane(self,
                      promissor: crd.Vector,
                      acceptor: crd.Vector,
                      aspect: Optional[int] = None) -> Optional[list[dict]]:
        """
        Returns mundane primary directions of the promissor
        to aspect points of the acceptor
        """
        # Find acceptor's aspect positions on the equator
        aspect_positions = self.aspect_positions_regio_mundane(
            acceptor, aspect)

        promissor_eqt = promissor.equatorial()
        promissor_rasc = promissor_eqt.rasc
        promissor_dec = promissor_eqt.dec

        return [
            dict(
                dist=promissor_rasc - self.__ra_from_oa_utp__(
                    promissor_dec, item['rasc']),
                aspect=item['aspect']
            ) for item in aspect_positions
        ]


# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
# Example of usage:
# ╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌╌
if __name__ == '__main__':
    # Initiate a celestial sphere object
    test_sphere = Sphere(
        datetime(2022, 10, 10, 23, 20),
        time_zone=4,
        geo_lon=44 + 46/60,
        geo_lat=56 + 43/60,
    )

    test_directions = Directions(test_sphere)

    # Set the 1st point
    proms = test_sphere.set_ecliptical(150, 0)

    # Set the 2nd point
    accpt = test_sphere.set_ecliptical(140, 0)

    print(test_directions.placidus_mundane(proms, accpt, aspect=None))
    print(test_directions.placidus_zodiac(proms, accpt, aspect=0))
    print(test_directions.placidus_zodiac(
        proms, accpt, aspect=0, field_plane_lat=10))
    print(test_directions.regio_mundane(proms, accpt, aspect=0))
