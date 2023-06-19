import React, { Component } from 'react';
import Map from 'ol/Map';
import View from 'ol/View';
import { Tile as TileLayer, Vector as VectorLyr, Image } from 'ol/layer';
import { OSM, Vector as VectorSrc, ImageCanvas } from 'ol/source';
import Feature from 'ol/Feature';
import { Point, Polygon } from 'ol/geom';
import { Style, Circle, Fill } from 'ol/style'

import './custom.css'
import lyrData from './mock/layer1.json';
import kriging from './kriging.js'

export default class App extends Component {
    componentDidMount() {
        const vctSource = new VectorSrc();
        const vctLayer = new VectorLyr({ source: vctSource });

        const map = new Map({
            target: "map",
            view: new View({
                center: [116.40, 39.90],
                projection: 'EPSG:4326',
                zoom: 12
            }),
            layers: [
                new TileLayer({
                    source: new OSM({
                        "url": "http://tile2.opencyclemap.org/transport/{z}/{x}/{y}.png"
                    })
                }),
                vctLayer
            ]
        })

        const { features, coord, colors } = lyrData;
        features.forEach((f) => {
            const feature = new Feature({
                geometry: new Point([f.geometry.x, f.geometry.y]),
                value: f.attributes.z,
            })
            feature.setStyle(
                new Style({
                    image: new Circle({
                        radius: 6,
                        fill: new Fill({
                            color: '#00F',
                        }),
                    }),
                }),
            )

            vctSource.addFeature(feature);
        });

        setTimeout(async () => {
            const points = features.map(f => ({ ...f.geometry, value: f.attributes.z }));

            const valueList = points.map(p => p.value).sort((a, b) => a - b);
            const zlim = [valueList[0], valueList[valueList.length - 1]];

            const polygons = [({ points: coord.map(([x, y]) => ({ x, y })) })];

            const exPolygon = new Polygon([coord]);
            const ex = exPolygon.getExtent();
            const width = (ex[2] - ex[0]) / 200;

            const response = await fetch('kriging', {
                method: 'post',
                headers: {
                    'Content-Type': 'application/json'
                },
                body: JSON.stringify({
                    points, sigma2: 0, alpha: 256, polygons, width
                })
            });
            const grid = await response.json();

            const canvasLayer = new Image({
                source: new ImageCanvas({
                    canvasFunction: (extent, resolution, pixelRatio, size, projection) => {
                        let canvas = document.createElement('canvas')
                        canvas.width = size[0]
                        canvas.height = size[1]
                        canvas.style.display = 'block'
                        canvas.getContext('2d').globalAlpha = 0.75

                        kriging.plot(canvas, {
                            ...grid.grid,
                            length: grid.grid.length,
                            xlim: grid.xlim,
                            ylim: grid.ylim,
                            zlim,
                            width
                        }, [extent[0], extent[2]], [extent[1], extent[3]], colors);


                        //  !!! code below is also valid ,but i wonder where to put 

                        //// plot 
                        //var ctx = canvas.getContext("2d");
                        //ctx.clearRect(0, 0, canvas.width, canvas.height);

                        //const xlim = [extent[0], extent[2]]
                        //const ylim = [extent[1], extent[3]]
                        //

                        //// Starting boundaries
                        //var range = [xlim[1] - xlim[0], ylim[1] - ylim[0], zlim[1] - zlim[0]];
                        //var i, j, x, y, z;

                        //var n = grid.grid.length;
                        //var m = grid.grid[0].length;
                        //var wx = Math.ceil(width * canvas.width / (xlim[1] - xlim[0]));
                        //var wy = Math.ceil(width * canvas.height / (ylim[1] - ylim[0]));
                        //for (i = 0; i < n; i++)
                        //    for (j = 0; j < m; j++) {
                        //        if (grid.grid[i][j] == undefined) continue;

                        //        x = canvas.width * (i * width + grid.xlim[0] - xlim[0]) / range[0];
                        //        y = canvas.height * (1 - (j * width + grid.ylim[0] - ylim[0]) / range[1]);
                        //        z = (grid.grid[i][j] - zlim[0]) / range[2];

                        //        if (z < 0.0) z = 0.0;
                        //        if (z > 1.0) z = 1.0;

                        //        ctx.fillStyle = colors[Math.floor((colors.length - 1) * z)];
                        //        ctx.fillRect(Math.round(x - wx / 2), Math.round(y - wy / 2), wx, wy);
                        //    }

                        return canvas
                    },
                    projection: 'EPSG:4326',
                })
            });
            map.addLayer(canvasLayer);
        }, 0);
    }

    render() {
        return (
            <div id="map" />
        );
    }
}
