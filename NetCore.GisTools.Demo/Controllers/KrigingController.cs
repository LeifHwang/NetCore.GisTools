using Microsoft.AspNetCore.Mvc;

using NetCore.GisTools.Model;

using System;
using System.Diagnostics;
using System.Threading.Tasks;

namespace NetCore.GisTools.Demo.Controllers
{
    [ApiController]
    [Route("[controller]")]
    public class KrigingController : ControllerBase
    {
        [HttpPost]
        public async Task<IActionResult> Post([FromBody] KrigingInputModel model)
        {
            Console.WriteLine();
            Console.WriteLine($"kriging @{DateTime.Now.TimeOfDay}");

            var stopwatch = new Stopwatch();

            // train
            stopwatch.Start();
            var variogram = await Kriging.TrainAsync(model, 0, 100);

            stopwatch.Stop();
            Console.WriteLine($"train cost: {stopwatch.Elapsed}");

            //grid
            stopwatch.Restart();
            var grid = Kriging.Grid(model, variogram, (model.MaxY - model.MinY) / 500);

            stopwatch.Stop();
            Console.WriteLine($"grid cost: {stopwatch.Elapsed}");

            return Ok(grid);
        }
    }
}
