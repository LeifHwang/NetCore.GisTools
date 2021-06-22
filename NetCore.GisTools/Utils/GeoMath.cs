using System;

namespace NetCore.GisTools.Utils
{
    public static class GeoMath
    {
        public static double Distance(double x1, double y1, double x2, double y2)
        {
            var x = x1 - x2;
            var y = y1 - y2;
            return Math.Sqrt((x * x) + (y * y));
        }
    }
}
