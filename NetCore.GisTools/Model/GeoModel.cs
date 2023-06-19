using System;
using System.Diagnostics;

namespace NetCore.GisTools.Model
{
    [DebuggerDisplay("({Longitude},{Latitude})")]
    public class LocationPoint
    {
        public double Longitude { get; set; }
        public double Latitude { get; set; }

        public LocationPoint() { }
        public LocationPoint(double x, double y)
        {
            Longitude = x;
            Latitude = y;
        }

        public double Distance<T2>(GeoData<T2> other)
        {
            return Math.Pow(Math.Pow(Longitude - other.Longitude, 2) +
                        Math.Pow(Latitude - other.Latitude, 2), 0.5);
        }
    }

    public class GeoPolygon
    {
        public GeoPolygon(LocationPoint[] points)
        {
            Points = points;
        }

        public LocationPoint[] Points { get; }

        public bool Pip(double x, double y) => Pip(new LocationPoint(x, y));
        public bool Pip(LocationPoint testPoitn)
        {
            for (int i = 0, j = Points.Length - 1; i < Points.Length; j = i++)
            {
                var current = Points[i];
                var previous = Points[j];

                if (((current.Latitude > testPoitn.Latitude) != (previous.Latitude > testPoitn.Latitude)) &&
                    (testPoitn.Longitude < ((previous.Longitude - current.Longitude) * (testPoitn.Latitude - current.Latitude) / (previous.Latitude - current.Latitude)) + current.Longitude))
                {
                    return true;
                }
            }
            return false;
        }
    }
}
