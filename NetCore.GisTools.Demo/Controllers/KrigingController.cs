using Microsoft.AspNetCore.Mvc;

using NetCore.GisTools.Demo.Model;
using NetCore.GisTools.Model;

using System.Linq;

namespace NetCore.GisTools.Demo.Controllers
{
    [ApiController]
    [Route("[controller]")]
    public class KrigingController : ControllerBase
    {
        [HttpPost]
        public object Demo([FromBody] KrigingInput data)
        {
            var variogram = Kriging.Train(
                data.points.Select(i => new GeoDoubleData(i.x, i.y, i.value)).ToArray(),
                data.sigma2,
                data.alpha);

            var polygons = data.polygons
                .Select(py => new GeoPolygon(py.points.Select(p => new GeoPoint(p.x, p.y)).ToArray())).ToArray();
            var gird = Kriging.Grid(polygons, variogram, data.width);

            var obj = gird.ForPlot(variogram);

            return Ok(obj);
        }

        
    }
}
