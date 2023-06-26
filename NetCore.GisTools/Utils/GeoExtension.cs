using NetCore.GisTools.Model;

using System;

namespace NetCore.GisTools.Utils
{
    static class GeoExtension
    {
        public static double Distance(this GeoPoint from, GeoPoint to) => 
            Math.Pow(Math.Pow(from.Longitude - to.Longitude, 2) + 
                Math.Pow(from.Latitude - to.Latitude, 2), 0.5);
    }
}
