# NetCore.GisTools
Tools for gis by the .net core, includes [kriging](#kriging) ...

### <a id="kriging"></a> `Kriging`
  > Refrence to [kriging.js](https://github.com/oeo4b/kriging.js)

Main code for usage:
- `var variogram = Kriging.Train(geoDataArray, sigma2 = 0, alpha = 100, model = KrigingTrainModel.Exponential);`
- `var gird = Kriging.Grid(polygons, variogram, gridWidth);`
- if plot by kriging.js, then call `var obj = gird.ForPlot(variogram);` to converter the grid to ParamreterType for `kriging.plot(canvas, grid, xlim, ylim, colors);`

  *Or u can also call `var list = Kriging.Plot(variogram, grid, options);` to get all pixel points that can be drawn on canvas*
