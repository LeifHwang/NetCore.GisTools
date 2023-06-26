using System;
using System.Diagnostics;

namespace NetCore.GisTools.Model
{
    [DebuggerDisplay("({Longitude},{Latitude})")]
    public class GeoPoint
    {
        public double Longitude { get; }
        public double Latitude { get; }

        public GeoPoint(double longitude, double latitude)
        {
            if (longitude < -180 || longitude > 180)
                throw new ArgumentOutOfRangeException(nameof(longitude), "invalid longitude");
            if (latitude < -90 || latitude > 90)
                throw new ArgumentOutOfRangeException(nameof(latitude), "invalid latitude");

            Longitude = longitude;
            Latitude = latitude;
        }
    }

    public class GeoPolygon
    {
        public GeoPoint[] Points { get; }

        public GeoPolygon(GeoPoint[] points)
        {
            if (points.Length < 3)
                throw new ArgumentOutOfRangeException(nameof(points), "At least three points are required!");

            Points = points;
        }


        public GeoExtend GetExtend()
        {
            double xMin = Points[0].Longitude, xMax = Points[0].Longitude;
            double yMin = Points[0].Latitude, yMax = Points[0].Latitude;

            for (int i = 1; i < Points.Length; i++)
            {
                var point = Points[i];

                xMin = Math.Min(xMin, point.Longitude);
                xMax = Math.Max(xMax, point.Longitude);

                yMin = Math.Min(yMin, point.Latitude);
                yMax = Math.Max(yMax, point.Latitude);
            }

            return new GeoExtend(xMin, xMax, yMin, yMax);
        }

        public bool Pip(double x, double y) => Pip(new GeoPoint(x, y));
        public bool Pip(GeoPoint testPoitn)
        {
            var result = false;
            for (int i = 0, j = Points.Length - 1; i < Points.Length; j = i++)
            {
                var x1 = Points[i].Longitude;
                var y1 = Points[i].Latitude;
                var x2 = Points[j].Longitude;
                var y2 = Points[j].Latitude;

                if (((y1 > testPoitn.Latitude) != (y2 > testPoitn.Latitude)) &&
                    (testPoitn.Longitude < (x2 - x1) * (testPoitn.Latitude - y1) / (y2 - y1) + x1))
                {
                    result = !result;
                }

            }
            return result;
        }
    }

    [DebuggerDisplay("{LngMin},{LngMax};{LatMin},{LatMax}")]
    public class GeoExtend
    {
        public double LngMin { get; private set; }
        public double LngMax { get; private set; }

        public double LatMin { get; private set; }
        public double LatMax { get; private set; }

        public double Width => LngMax - LngMin;
        public double Height => LatMax - LatMin;

        public GeoExtend(double xMin, double xMax, double yMin, double yMax)
        {
            LngMin = xMin;
            LngMax = xMax;

            LatMin = yMin;
            LatMax = yMax;
        }

        /// <summary>
        /// merge other rectangle
        /// </summary>
        /// <remarks>will expand current extend when the other rec intersect or outside</remarks>
        public void Union(GeoExtend other)
        {
            LngMin = Math.Min(LngMin, other.LngMin);
            LngMax = Math.Max(LngMax, other.LngMax);

            LatMin = Math.Min(LatMin, other.LatMin);
            LatMax = Math.Max(LatMax, other.LatMax);
        }
    }
}
