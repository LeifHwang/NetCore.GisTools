using NetCore.GisTools.Model;

namespace NetCore.GisTools.Utils
{
    public static class GeoDataExtensions
    {
        public static double CalcDistance<T>(this GeoData<T> data, double x, double y)
        {
            return GeoMath.Distance(data.Longitude, x, data.Latitude, y);
        }

        public static double DistanceBetween<T>(this GeoData<T> from, GeoData<T> to)
        {
            return from.CalcDistance(to.Longitude, to.Latitude);
        }
    }
}
