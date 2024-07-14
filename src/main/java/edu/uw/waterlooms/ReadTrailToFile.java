package edu.uw.waterlooms;

import edu.uw.waterlooms.entity.IsolationWindow;
import edu.uw.waterlooms.entity.XIC;
import edu.uw.waterlooms.match.AminoAcid;
import me.tongfei.progressbar.ProgressBar;
import org.json.JSONWriter;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;

public class ReadTrailToFile {
    public class MSTwoWindow{
        public double mzHigh;
        public double mzLow;
    }
    public List<IsolationWindow> windows = new ArrayList<>();

    public static String html_head = "<!--\n" +
            "    THIS EXAMPLE WAS DOWNLOADED FROM https://echarts.apache.org/examples/zh/editor.html?c=grid-multiple\n" +
            "-->\n" +
            "<!DOCTYPE html>\n" +
            "<html style='height: 100%'>\n" +
            "\n" +
            "<head>\n" +
            "    <meta charset='utf-8'>\n" +
            "</head>\n" +
            "\n" +
            "<body style='height: 100%; margin: 0'>\n" +
            "    <div  style='height: 5%'>\n" +
            "        Enter peptide: <input type='text' id='myText' name='txt'  onchange='myFunction(this.value)'>\n" +
            "\n" +
            "    </div>\n" +
            "    <div id='container' style='height: 95%'>\n" +
            "\n" +
            "    </div>\n" +
            "    <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/dist/echarts.min.js'></script>\n" +
            "    <script type='text/javascript' src='";
    public static String html_JS = "'></script>\n" +
            "\n" +
            "    <!-- Uncomment this line if you want to dataTool extension\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/dist/extension/dataTool.min.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment this line if you want to use gl extension\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts-gl@2/dist/echarts-gl.min.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment this line if you want to echarts-stat extension\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts-stat@latest/dist/ecStat.min.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment this line if you want to use map\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/map/js/china.js'></script>\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@5.3.0/map/js/world.js'></script>\n" +
            "        -->\n" +
            "    <!-- Uncomment these two lines if you want to use bmap extension\n" +
            "        <script type='text/javascript' src='https://api.map.baidu.com/api?v=2.0&ak=<Your Key Here>'></script>\n" +
            "        <script type='text/javascript' src='https://cdn.jsdelivr.net/npm/echarts@{{version}}/dist/extension/bmap.min.js'></script>\n" +
            "        -->\n" +
            "    <script type='text/javascript'>\n" +
            "    var dom = document.getElementById('container');\n" +
            "    var myChart = echarts.init(dom);\n" +
            "    var app = {};\n" +
            "\n" +
            "    var option;\n" +
            "    var differTwoIons = 0.0;\n" +
            "    var TimeLinecurrentIndex = 0;\n" +
            "\n" +
            "    var ProtonMass = 1.00727647;\n" +
            "    var HMass = 1.007825035;\n" +
            "    var H2OMass = 18.0105647;\n" +
            "    const aminoAcids = new Map();\n" +
            "\n" +
            "    aminoAcids.set('A', 71.037114);\n" +
            "    aminoAcids.set('R', 156.101111);\n" +
            "    aminoAcids.set('N', 114.042927);\n" +
            "    aminoAcids.set('D', 115.026943);\n" +
            "    aminoAcids.set('C', 160.030649);//103.009185);//160.03065 add fixed ptm 二琉键 57.021464\n" +
            "    aminoAcids.set('Q', 128.058578);\n" +
            "    aminoAcids.set('E', 129.042593);\n" +
            "    aminoAcids.set('G', 57.021464);\n" +
            "    aminoAcids.set('H', 137.058912);\n" +
            "    aminoAcids.set('I', 113.084064);\n" +
            "    aminoAcids.set('L', 113.084064);\n" +
            "    aminoAcids.set('K', 128.094963);\n" +
            "    aminoAcids.set('M', 131.040485);\n" +
            "    aminoAcids.set('F', 147.068414);\n" +
            "    aminoAcids.set('P', 97.052764);\n" +
            "    aminoAcids.set('S', 87.032028);\n" +
            "    aminoAcids.set('T', 101.047679);\n" +
            "    aminoAcids.set('W', 186.079313);\n" +
            "    aminoAcids.set('Y', 163.063332);\n" +
            "    aminoAcids.set('V', 99.068414);\n" +
            "\n" +
            "    var peptide = '";
    public static String html_peptide = "';\n" +
            "    var len = peptide.length;\n" +
            "\n" +
            "    var bionsMasses = new Array(len);\n" +
            "    var b2ChargeionsMasses = new Array(len);\n" +
            "    var yionsMasses = new Array(len);\n" +
            "    var y2ChargeionsMasses = new Array(len);\n" +
            "\n" +
            "    var PeptideMasssum = 18.0105647;\n" +
            "    var partial_y = H2OMass + ProtonMass ;\n" +
            "    var partial_b = ProtonMass;\n" +
            "\n" +
            "    // var b1Intensity = [0,0.49706772,0.149048626,0.31874302,0.125484154,0.114351839,0.0498381,0.050654951,0.011937074,0.013909766,0.000594605,0,0,0,0.001573052,0];\n" +
            "    // var b2Intensity = [0,0,0,0.000948123,0,0,0,0,0,0,0,0,0,0,0,0];\n" +
            "    // var y1Intensity = [0.061157208,0.559807897,0.075900644,0.431469232,0.204922244,0.174942195,0.335118979,0.274580389,1,0.324403346,0.836681426,0.451821238,0.174080625,0.123319842,0.102458522,0];\n" +
            "    // var y2Intensity = [0,0,0,0,0,0,0,0,0,0,0,0.002780302,0.021060117,0,0.037728656,0];\n" +
            "\n" +
            "    var Intensity= [[2.47084881e-14,0.00000000e+00,5.75685488e-18,2.22559904e-09]\n" +
            ",[2.34642539e-05,1.30377115e-11,3.54840376e-15,7.91612183e-05]\n" +
            ",[1.63372070e-01,1.16119480e-07,1.40902711e-13,3.55197638e-02]\n" +
            ",[1.46155819e-01,1.82013096e-07,1.01242463e-12,8.26046690e-02]\n" +
            ",[1.63832679e-01,2.36935787e-08,5.20611852e-12,1.52841330e-01]\n" +
            ",[2.23015890e-01,5.25931410e-10,5.02555159e-12,2.03018636e-01]\n" +
            ",[1.98946223e-01,7.41375272e-09,4.16820076e-05,9.98054683e-01]\n" +
            ",[2.65103090e-03,1.24728017e-12,1.91180340e-12,2.96788141e-02]\n" +
            ",[6.02396764e-03,2.35211054e-14,1.27342332e-08,4.42096218e-02]\n" +
            ",[1.07828579e-04,2.36312353e-13,1.22322710e-02,8.75109285e-02]\n" +
            ",[4.93524969e-12,2.72961788e-14,1.78765487e-02,2.78820805e-02]\n" +
            ",[1.44625312e-09,8.02632537e-13,6.04677871e-02,6.05705492e-02]\n" +
            ",[2.72290409e-03,2.79258017e-09,3.90534848e-01,1.67324454e-01]\n" +
            ",[0.00000000e+00,8.76269435e-13,3.06502972e-02,1.24782976e-03]\n" +
            ",[4.71512423e-13,4.35996031e-11,6.94957227e-02,1.42859966e-02]\n" +
            ",[1.22912357e-11,5.01179273e-15,1.25387773e-01,7.56352162e-03]\n" +
            ",[0.00000000e+00,0.00000000e+00,2.77467668e-02,1.32186706e-09]\n" +
            ",[0.00000000e+00,0.00000000e+00,5.48798479e-02,1.61414815e-09]\n" +
            ",[0.00000000e+00,0.00000000e+00,1.87662035e-01,1.95182941e-03]\n" +
            ",[0.00000000e+00,0.00000000e+00,9.88203567e-03,3.77472457e-18]\n" +
            ",[0.00000000e+00,0.00000000e+00,8.05829372e-03,8.43143573e-22]\n" +
            ",[0.00000000e+00,0.00000000e+00,2.76528709e-02,1.35532651e-26]\n" +
            ",[0.00000000e+00,0.00000000e+00,3.86985317e-02,1.84890576e-24]\n" +
            ",[0.00000000e+00,0.00000000e+00,3.66997831e-02,1.09663568e-24]\n" +
            ",[0.00000000e+00,0.00000000e+00,1.31967507e-04,5.31150053e-29]\n" +
            ",[0.00000000e+00,0.00000000e+00,1.10709511e-01,1.55309446e-10]\n" +
            ",[0.00000000e+00,0.00000000e+00,4.19339376e-05,6.20807218e-19]];\n" +
            "\n" +
            "    var intensitybz = false;\n" +
            "\n" +
            "    for (var i = 0; i < peptide.length-1; i++) {\n" +
            "\n" +
            "          \n" +
            "          PeptideMasssum += aminoAcids.get(peptide.charAt(i));\n" +
            "          partial_y += aminoAcids.get(peptide.charAt(len - i -1));\n" +
            "          partial_b += aminoAcids.get(peptide.charAt(i));\n" +
            "          bionsMasses[i] = new Array(3);\n" +
            "          bionsMasses[i][0] = partial_b; \n" +
            "          if (intensitybz)\n" +
            "            bionsMasses[i][1] = 100 * Intensity[i][0]; \n" +
            "            \n" +
            "          else \n" +
            "            bionsMasses[i][1] = 100;\n" +
            "          bionsMasses[i][2] = i+1; \n" +
            "\n" +
            "          b2ChargeionsMasses[i] = new Array(3);\n" +
            "          b2ChargeionsMasses[i][0] = (partial_b + ProtonMass)/2; \n" +
            "          if (intensitybz)\n" +
            "            b2ChargeionsMasses[i][1] = 100 * Intensity[i][1]; \n" +
            "          else\n" +
            "            b2ChargeionsMasses[i][1] = 90\n" +
            "          b2ChargeionsMasses[i][2] = i+1; \n" +
            "\n" +
            "          yionsMasses[i] = new Array(3);\n" +
            "          yionsMasses[i][0] = partial_y; \n" +
            "          if (intensitybz)\n" +
            "            yionsMasses[i][1] = 100 * Intensity[i][2]; \n" +
            "          else\n" +
            "            yionsMasses[i][1] = 100\n" +
            "          yionsMasses[i][2] = i+1; \n" +
            "\n" +
            "          y2ChargeionsMasses[i] = new Array(3);\n" +
            "          y2ChargeionsMasses[i][0] = (partial_y + ProtonMass)/2;\n" +
            "          if (intensitybz)\n" +
            "            y2ChargeionsMasses[i][1] = 100 * Intensity[i][3];\n" +
            "          else\n" +
            "            y2ChargeionsMasses[i][1] = 90\n" +
            "          y2ChargeionsMasses[i][2] = i+1; \n" +
            "\n" +
            "    }\n" +
            "\n" +
            "\n" +
            "    // prettier-ignore\n" +
            "    // let timeData = [\n" +
            "    //     '2009/10/18 8:00'\n" +
            "    // ];\n" +
            "    // timeData = timeData.map(function (str) {\n" +
            "    //   return str.replace('2009/', '');\n" +
            "    // });\n" +
            "    document.getElementById(\"myText\").value = peptide;\n" +
            "    option = {\n" +
            "     baseOption: {\n" +
            "\n" +
            "        timeline: {\n" +
            "          position: 'top',\n" +
            "\n" +
            "          axisType: 'category',\n" +
            "          // realtime: false,\n" +
            "          // loop: false,\n" +
            "          // autoPlay: true,\n" +
            "          // currentIndex: 2,\n" +
            "          // playInterval: 2000,\n" +
            "          // controlStyle: {\n" +
            "          //     position: 'left'\n" +
            "          // },\n" +
            "          data: ms2RTRange\n" +
            "          \n" +
            "        },\n" +
            "\n" +
            "        title: {\n" +
            "            text: 'BYIons vs ALL PEAKs_'+peptide+'_MS1RT:'+ms1RT+'_MS2RT:'+ms2RT+'_'+PeptideMasssum,\n" +
            "            left: 'center'\n" +
            "        },\n" +
            "        tooltip: {\n" +
            "            trigger: 'axis',\n" +
            "            axisPointer: {\n" +
            "                animation: true\n" +
            "            }\n" +
            "        },\n" +
            "        legend: {\n" +
            "            data: ['MS2Trail', 'B', 'B+', 'Y', 'Y+'],\n" +
            "            left: 10\n" +
            "        },\n" +
            "        toolbox: {\n" +
            "            feature: {\n" +
            "                dataZoom: {\n" +
            "                    yAxisIndex: 'none'\n" +
            "                },\n" +
            "                restore: {},\n" +
            "                saveAsImage: {}\n" +
            "            }\n" +
            "        },\n" +
            "        axisPointer: {\n" +
            "            link: [{\n" +
            "                xAxisIndex: 'all'\n" +
            "            }]\n" +
            "        },\n" +
            "        dataZoom: [{\n" +
            "                show: true,\n" +
            "                realtime: true,\n" +
            "                start: 0,\n" +
            "                end: 100,\n" +
            "                xAxisIndex: [0, 1]\n" +
            "            },\n" +
            "            {\n" +
            "                type: 'inside',\n" +
            "                realtime: true,\n" +
            "                start: 30,\n" +
            "                end: 70,\n" +
            "                xAxisIndex: [0, 1]\n" +
            "            }\n" +
            "        ],\n" +
            "        grid: [{\n" +
            "                left: 80,\n" +
            "                right: 50,\n" +
            "                top: '8%',\n" +
            "\n" +
            "                height: '40%'\n" +
            "            },\n" +
            "            {\n" +
            "                left: 80,\n" +
            "                right: 50,\n" +
            "                top: '48%',\n" +
            "                height: '40%'\n" +
            "            }\n" +
            "        ],\n" +
            "        xAxis: [{\n" +
            "                position: 'top',\n" +
            "                type: 'value',\n" +
            "                min: 0,\n" +
            "                max: partial_y+100,//2103.8067367,\n" +
            "                //   type: 'category',\n" +
            "                boundaryGap: false,\n" +
            "                axisPointer: {\n" +
            "                    label: {\n" +
            "                      formatter: function (params) {\n" +
            "                         if (typeof params.seriesData[0] !== 'undefined' && typeof params.seriesData[0].data !== 'undefined')\n" +
            "                        {\n" +
            "                            differTwoIons =  params.seriesData[0].data[0] ;\n" +
            "                        \n" +
            "\n" +
            "                            return (\n" +
            "                              'mass  ' +\n" +
            "                              params.value  +\n" +
            "                              (typeof params.seriesData[0].data !== 'undefined' ? ' NO.：' + params.seriesData[0].data[2] : '')\n" +
            "                            );\n" +
            "                        }\n" +
            "                      }\n" +
            "                    }\n" +
            "                },\n" +
            "                axisLine: { onZero: true } //,\n" +
            "                //   data: timeData\n" +
            "            },\n" +
            "            {\n" +
            "                gridIndex: 1,\n" +
            "                //     type: 'category',\n" +
            "                type: 'value',\n" +
            "                min: 0,\n" +
            "                max:partial_y+100,//2103.8067367,\n" +
            "                boundaryGap: false,\n" +
            "                axisPointer: {\n" +
            "                    label: {\n" +
            "                      formatter: function (params) {\n" +
            "\n" +
            "                        return (\n" +
            "                          'mass  ' +\n" +
            "                          params.value  +'\\n'+\n" +
            "                          'diff:'+(Math.abs(params.value-differTwoIons))+\n" +
            "                          (((Math.abs(params.value-differTwoIons)*1000000)/params.value)<10 ? ' Matched': 'NoMatch')\n" +
            "                        );\n" +
            "                      }\n" +
            "                    }\n" +
            "                },\n" +
            "                axisLine: { onZero: true } //,\n" +
            "\n" +
            "                //      data: timeData\n" +
            "            }\n" +
            "        ],\n" +
            "        yAxis: [{\n" +
            "                name: 'MS2Trail',\n" +
            "                type: 'value'\n" +
            "                // ,\n" +
            "                // max: 2000\n" +
            "            },\n" +
            "            {\n" +
            "                gridIndex: 1,\n" +
            "                name: 'BYIons',\n" +
            "                type: 'value',\n" +
            "                inverse: true\n" +
            "                 ,\n" +
            "                max: 110\n" +
            "            }\n" +
            "        ],\n" +
            "        series: [\n" +
            "\n" +
            "            {\n" +
            "                name: 'MS2Trail',\n" +
            "                type: 'bar',\n" +
            "                symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#0001ff'\n" +
            "                    }\n" +
            "                },\n" +
            "                // label: {\n" +
            "                //   normal: {\n" +
            "                //     show: true,\n" +
            "                //     color: '#323232',\n" +
            "                //     position: 'top'\n" +
            "                //   }\n" +
            "                // },\n" +
            "                // prettier-ignore\n" +
            "                data: ms2data[0]\n" +
            "            },\n" +
            "            {\n" +
            "                name: 'B',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                        // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'blue',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#49a1ff'\n" +
            "                    }\n" +
            "                },\n" +
            "                // label: {\n" +
            "                //   normal: {\n" +
            "                //     show: true,\n" +
            "                //     color: '#323232',\n" +
            "                //     position: 'top'\n" +
            "                //   }\n" +
            "                // },\n" +
            "                // prettier-ignore\n" +
            "                data:bionsMasses\n" +
            "            },\n" +
            "            {\n" +
            "                name: 'B+',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#ff00ff'\n" +
            "                    }\n" +
            "                },\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                        // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'red',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                // prettier-ignore\n" +
            "                data: b2ChargeionsMasses\n" +
            "            },\n" +
            "\n" +
            "            {\n" +
            "                name: 'Y',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#4ff0a1'\n" +
            "                    }\n" +
            "                },\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                        // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'red',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                // prettier-ignore\n" +
            "\n" +
            "                data:yionsMasses\n" +
            "            },\n" +
            "            {\n" +
            "                name: 'Y+',\n" +
            "                type: 'bar',\n" +
            "                xAxisIndex: 1,\n" +
            "                yAxisIndex: 1,\n" +
            "                // symbolSize: 8,\n" +
            "                barWidth: '1px',\n" +
            "\n" +
            "                itemStyle: {\n" +
            "                    normal: {\n" +
            "                        color: '#ff8f00'\n" +
            "                    }\n" +
            "                },\n" +
            "                label: {\n" +
            "                    normal: {\n" +
            "                        show: true,\n" +
            "                        color: '#323232',\n" +
            "                        position: 'bottom',\n" +
            "                       // formatter: '{a|{a}  \\n{@[2]}}', //\\n{c}',\n" +
            "                        //\n" +
            "                        //  function(e){\n" +
            "                        //   let data=e.data;\n" +
            "                        //   return '{a}${data.value}:+{a|这段文本采用样式a}';\n" +
            "                        // },\n" +
            "\n" +
            "                        formatter: function(d) {\n" +
            "                            // console.log(d);\n" +
            "                                return formatlab(d);\n" +
            "                              },\n" +
            "\n" +
            "\n" +
            "\n" +
            "                        rich: {\n" +
            "                            a: {\n" +
            "                                align: 'center',\n" +
            "                                color: 'blue',\n" +
            "                                lineHeight: 10\n" +
            "                            },\n" +
            "                            b: {\n" +
            "                                backgroundColor: {\n" +
            "                                    image: 'xxx/xxx.jpg'\n" +
            "                                },\n" +
            "                                height: 40\n" +
            "                            },\n" +
            "                            x: {\n" +
            "                                // fontSize: 18,\n" +
            "                                // fontFamily: 'Microsoft YaHei',\n" +
            "                                align: 'center',\n" +
            "\n" +
            "                                borderColor: '#449933'\n" +
            "                                // ,\n" +
            "                                // borderWidth: 1,\n" +
            "\n" +
            "                                // borderRadius: 1\n" +
            "                            },\n" +
            "\n" +
            "                        }\n" +
            "                    },\n" +
            "\n" +
            "                },\n" +
            "                // prettier-ignore\n" +
            "                // \n" +
            "                data: y2ChargeionsMasses\n" +
            "\n" +
            "            }\n" +
            "        ]\n" +
            "        },\n" +
            "        options: optionsText\n" +
            "    };\n" +
            "\n" +
            "    if (option && typeof option === 'object') {\n" +
            "        myChart.setOption(option);\n" +
            "    }\n" +
            "function isMatched(dataMass)\n" +
            "    {\n" +
            "        //  console.log(option.options[TimeLinecurrentIndex]);\n" +
            "        // console.log(option.baseOption.timeline);  \n" +
            "        var dinteralDiffThreshold = dataMass * 10/1000000;\n" +
            "        var bMatched  = iterativeFunction(option.options[TimeLinecurrentIndex].series[0].data,dataMass,dinteralDiffThreshold);\n" +
            "\n" +
            "        return bMatched;\n" +
            "\n" +
            "    }\n" +
            "    function iterativeFunction (arr, x, interalvalue) {\n" +
            "  \n" +
            "        var start=0, end=arr.length-1;\n" +
            "\n" +
            "             \n" +
            "        // Iterate while start not meets end\n" +
            "        while (start<=end){\n" +
            "             \n" +
            "\n" +
            "            // Find the mid index\n" +
            "            var mid=Math.floor((start + end)/2);\n" +
            "      \n" +
            "            // serial data is 3 array, If element is present at mid, return True\n" +
            "            // if (arr[0][mid]===x) return mid;\n" +
            "            if (Math.abs(arr[mid][0]-x)<interalvalue) return mid;\n" +
            "            // Else look in left or right half accordingly\n" +
            "            else if (arr[mid][0] < x)\n" +
            "                 start = mid + 1;\n" +
            "            else\n" +
            "                 end = mid - 1;\n" +
            "        }\n" +
            "      \n" +
            "        return -start;\n" +
            "    }\n" +
            "    function formatlab(d)\n" +
            "    {\n" +
            "        var match_h2o = isMatched(d.data[0]-18.0105647);\n" +
            "        var match_h = isMatched(d.data[0]-1.007825035);\n" +
            "        var match_nh3 = isMatched(d.data[0]-17.026549105);\n" +
            "\n" +
            "        return '{a|'+d.seriesName +' \\n'+d.data[2]+'}{x|'\n" +
            "        +((isMatched(d.data[0])>0)?'\\nMatch':'')\n" +
            "        // +((match_h2o>0)?('\\nM-H2O pos:'+(match_h2o+1)):'')\n" +
            "        // +((match_h>0)?('\\nM-H pos:'+(match_h+1)):'')\n" +
            "        // +((match_nh3>0)?('\\nM-NH3 pos:'+(match_nh3+1)):'')\n" +
            "        +'}';\n" +
            "    }    \n" +
            "\n" +
            "    //自动设置\n" +
            "// setInterval(function () {\n" +
            "//   console.log(option.dataset.source);\n" +
            "//   option.dataset.source[1][1] = +(Math.random() * 100).toFixed(2);\n" +
            "//   option.dataset.source[2][1] = +(Math.random() * 100).toFixed(2);\n" +
            "//   option.dataset.source[3][1] = +(Math.random() * 100).toFixed(2);\n" +
            "//   if (option && typeof option === 'object') {\n" +
            "//     myChart.setOption(option);\n" +
            "// }\n" +
            "//   // myChart.setOption({\n" +
            "//   //   series: [\n" +
            "//   //     {\n" +
            "//   //       data: option.dataset\n" +
            "//   //     }\n" +
            "//   //   ]\n" +
            "//   // });\n" +
            "// }, 3000);\n" +
            "function myFunction(val) {\n" +
            "  alert(\"The peptide value is: \" + val);\n" +
            "  // for (var i = 0; i < val.length; i++) {\n" +
            "\n" +
            "          \n" +
            "  //         partial_y += aminoAcids.get(peptide.charAt(len - i -1));\n" +
            "  //         partial_b += aminoAcids.get(peptide.charAt(i));\n" +
            "          \n" +
            "  //         option.options[TimeLinecurrentIndex].series[1].data[i][0] = partial_b; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[1].data[i][1] = 500; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[1].data[i][2] = i; \n" +
            "\n" +
            "          \n" +
            "  //         option.options[TimeLinecurrentIndex].series[2].data[i][0] = (partial_b + ProtonMass)/2; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[2].data[i][1] = 400; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[2].data[i][2] = i; \n" +
            "\n" +
            "         \n" +
            "  //         option.options[TimeLinecurrentIndex].series[3].data[i][0] = partial_y; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[3].data[i][1] = 500; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[3].data[i][2] = i+1; \n" +
            "\n" +
            "  //         option.options[TimeLinecurrentIndex].series[4].data[i][0] = (partial_y + ProtonMass)/2;\n" +
            "  //         option.options[TimeLinecurrentIndex].series[4].data[i][1] = 500; \n" +
            "  //         option.options[TimeLinecurrentIndex].series[4].data[i][2] = i+1; \n" +
            "\n" +
            "  //   }\n" +
            "  //   \n" +
            "  //   \n" +
            "    var lenChange = val.length;\n" +
            "    var bionsMassesChange = new Array(lenChange);\n" +
            "    var b2ChargeionsMassesChange = new Array(lenChange);\n" +
            "    var yionsMassesChange = new Array(lenChange);\n" +
            "    var y2ChargeionsMassesChange = new Array(lenChange);\n" +
            "\n" +
            "   \n" +
            "    var partial_yChange = H2OMass + ProtonMass ;\n" +
            "    var partial_bChange = ProtonMass;\n" +
            "\n" +
            "    for (var i = 0; i < lenChange; i++) {\n" +
            "\n" +
            "          \n" +
            "          partial_yChange += aminoAcids.get(val.charAt(lenChange - i -1));\n" +
            "          partial_bChange += aminoAcids.get(val.charAt(i));\n" +
            "          bionsMassesChange[i] = new Array(3);\n" +
            "          bionsMassesChange[i][0] = partial_bChange; \n" +
            "          bionsMassesChange[i][1] = 500; \n" +
            "          bionsMassesChange[i][2] = i+1; \n" +
            "\n" +
            "          b2ChargeionsMassesChange[i] = new Array(3);\n" +
            "          b2ChargeionsMassesChange[i][0] = (partial_bChange + ProtonMass)/2; \n" +
            "          b2ChargeionsMassesChange[i][1] = 400; \n" +
            "          b2ChargeionsMassesChange[i][2] = i+1; \n" +
            "\n" +
            "          yionsMassesChange[i] = new Array(3);\n" +
            "          yionsMassesChange[i][0] = partial_yChange; \n" +
            "          yionsMassesChange[i][1] = 500; \n" +
            "          yionsMassesChange[i][2] = i+1; \n" +
            "\n" +
            "          y2ChargeionsMassesChange[i] = new Array(3);\n" +
            "          y2ChargeionsMassesChange[i][0] = (partial_yChange + ProtonMass)/2;\n" +
            "          y2ChargeionsMassesChange[i][1] = 400; \n" +
            "          y2ChargeionsMassesChange[i][2] = i+1; \n" +
            "\n" +
            "    }\n" +
            "    if (option && typeof option === 'object') {\n" +
            "        option.options[TimeLinecurrentIndex].series[1].data = bionsMassesChange;\n" +
            "        option.options[TimeLinecurrentIndex].series[2].data = b2ChargeionsMassesChange;\n" +
            "        option.options[TimeLinecurrentIndex].series[3].data = yionsMassesChange;\n" +
            "        option.options[TimeLinecurrentIndex].series[4].data = y2ChargeionsMassesChange;\n" +
            "       \n" +
            "        myChart.setOption(option);\n" +
            "    }\n" +
            "  \n" +
            "    \n" +
            "}\n" +
            "myChart.on('timelinechanged', function(params) {\n" +
            "    TimeLinecurrentIndex =  params.currentIndex;\n" +
            "        console.log(TimeLinecurrentIndex);\n" +
            "    if (option && typeof option === 'object') {\n" +
            "        myChart.setOption(option);\n" +
            "    }\n" +
            "   \n" +
            "});\n" +
            "</script>\n" +
            "\n" +
            "</body>\n" +
            "\n" +
            "</html>";



    public void OutputMS2TrailInfo(String dataDirectory,String mzxmlfilename,String psmFilename,String outputPath) throws IOException {


        long start = System.currentTimeMillis();


        ArrayList<String> arrPepStr = new ArrayList<>();
        ArrayList<Double> arrdmz = new ArrayList<>();
        ArrayList<Double> arrdrt = new ArrayList<>();

        try {
            getPSMinfoFromFile(dataDirectory + psmFilename,arrPepStr,arrdmz,arrdrt);
//            IsolationWindowCollection(dataDirectory + "isolationWindowRanges.out");
        } catch (IOException e) {
            e.printStackTrace();
        }


        GenerateMS2Trail gen = new GenerateMS2Trail();

        windows = gen.ReadfeatureDetectShouldReturnAnArrayListOfXICs(dataDirectory, mzxmlfilename);//directory and mzXML filename


        double thresholdRT = 15.0/60;//5s
        ProgressBar pb = new ProgressBar("Progress to write file", arrPepStr.size());
        pb.start();
        for(int iPep=0;iPep<arrPepStr.size();iPep++)
        {
            pb.step();

            double dmz=arrdmz.get(iPep);//458.2248;
            double ms1rtTime =arrdrt.get(iPep);// 13.4918;
            boolean bzGetPeakArea =false;


            double ms1rtTimeBegin = ms1rtTime- thresholdRT*1.5;
            double ms1rtTimeEnd = ms1rtTime+ thresholdRT*1.5;




            int iwindows = FindWindowWithMZ(dmz,1);


//            System.out.println(dmz+" corresponding window "+iwindows);

            String fileWriterMS2TrailName =  "ms2trial"+ms1rtTimeBegin+"_"+ms1rtTimeEnd+".js";

            //得到以时间排序的treemap和以mz排序的list
            Map<Double, List<XIC>> mapRTMStwoTrailSorted =  gen.mapWindowXICs.get(iwindows).stream().sorted(XIC::compareMzXIC)
                    .collect(Collectors.groupingBy(XIC::getRtAtMaxIntensity, TreeMap::new, Collectors.toList()));
            FileWriter fileWriterMS2Trail = null;
            FileWriter fileWriterhtml = null;

//            String tmpPath = "ms2predictHtml/decoy/";
            try {
                fileWriterMS2Trail = new FileWriter(dataDirectory+outputPath+ fileWriterMS2TrailName);
                fileWriterhtml =  new FileWriter(dataDirectory+outputPath+iPep+"_"+arrPepStr.get(iPep)+ ".html");
            } catch (IOException e) {
                e.printStackTrace();
            }
            BufferedWriter bufferedWriterTest  = new BufferedWriter(fileWriterMS2Trail);
            BufferedWriter bufferedWriterHtml  = new BufferedWriter(fileWriterhtml );

            try {
                bufferedWriterTest.write("var ms2data = new Array();\n");
            } catch (IOException e) {
                e.printStackTrace();
            }



            final String[] RT = {"['",""};
            RT[1]="var optionsText = [\n" +
                    "            {\n" +
                    "              title: { text: '";
            AtomicInteger i = new AtomicInteger();
            mapRTMStwoTrailSorted.forEach((x,y)->{
//                System.out.println(x+" corresponding y  "+y.size()+" "+y.get(0));

                try {
                    if(x>=ms1rtTimeBegin && x<=ms1rtTimeEnd)
                    {
                        RT[0] = RT[0] +x+"','";
                        RT[1] = RT[1] +x+"' },\n" +
                                " series: [\n" +
                                "{ data: ms2data["+i+"]}]\n" +
                                "},\n" +
                                "{\n" +
                                "              title: { text: '";


                        bufferedWriterTest.write(" ms2data["+i+"] = ");
                        double[][] ddData= new double[y.size()][3];

                        int iL =0;
                        for(XIC msTwoTrail : y)
                        {
                            ddData[iL][0]=msTwoTrail.getMzAtMaxIntensity();

                            if(bzGetPeakArea) {
                                ddData[iL][1]=msTwoTrail.getPeakArea();
                            } else {
                                ddData[iL][1]=msTwoTrail.getCurIntensity(x);
                            }
                            ddData[iL][2]=iL+1;
                            iL++;
                        }
                        bufferedWriterTest.write(JSONWriter.valueToString(ddData)+";\n");
                        i.getAndIncrement();



                    }

                } catch (IOException e) {
                    e.printStackTrace();
                }
            });
            try {
                bufferedWriterTest.write("var ms2RT = "+ms1rtTimeBegin+";\n");
                bufferedWriterTest.write("var ms1RT = "+ms1rtTimeEnd+";\n");
                bufferedWriterTest.write("var ms2RTRange ="+RT[0].substring(0,RT[0].length()-2)+"];\n");
                bufferedWriterTest.write(RT[1].substring(0,RT[1].length()-34)+"];\n");

                bufferedWriterTest.flush();
                bufferedWriterTest.close();
                fileWriterMS2Trail.close();

                bufferedWriterHtml.write(html_head);
                bufferedWriterHtml.write(fileWriterMS2TrailName);
                bufferedWriterHtml.write(html_JS);
                bufferedWriterHtml.write(arrPepStr.get(iPep));
                bufferedWriterHtml.write(html_peptide);
                bufferedWriterHtml.flush();
                bufferedWriterHtml.close();
                fileWriterhtml.close();


            } catch (IOException e) {
                e.printStackTrace();
            }
        }
        pb.stop();
    }

    private void getPSMinfoFromFile(String sfilename, ArrayList<String> arrPepStr, ArrayList<Double> arrdmz, ArrayList<Double> arrdrt) throws IOException {

        String dataFile = sfilename;
        FileReader freader;
        try {
            freader = new FileReader(dataFile);
        } catch (FileNotFoundException noFile) {
            throw new FileNotFoundException();
        }

        try (BufferedReader br = new BufferedReader(freader)) {
            // create new isolation window
            String line;
            br.readLine();//skip titile
            while ((line = br.readLine()) != null) {
                // if line empty, continue
                if (line.isEmpty()) {
                    continue;
                }
                // if see end, add cur to list, next
                String[] psmInfo = line.split("\\s+");
                arrPepStr.add(psmInfo[0]);
                arrdmz.add(Double.parseDouble(psmInfo[1]));
                arrdrt.add(Double.parseDouble(psmInfo[2]));


            }
            br.close();
            freader.close();
        }
    }

//
//    public void IsolationWindowCollection(String trailInfile) throws IOException {
//
//        String dataFile = trailInfile;
//        FileReader freader;
//        try {
//            freader = new FileReader(dataFile);
//        } catch (FileNotFoundException noFile) {
//            throw new FileNotFoundException();
//        }
//
//        try (BufferedReader br = new BufferedReader(freader)) {
//            // create new isolation window
//            MSTwoWindow cur = new MSTwoWindow();
//            String line;
//
//            int lineType = 0; // 0: local mz/max I; 1: intensities of I; 2: retention times of trail
//
//            while ((line = br.readLine()) != null) {
//                // if line empty, continue
//                if (line.isEmpty()) {
//                    continue;
//                }
//                // if see end, add cur to list, next
//                if (line.contains("END")) {
//                    this.windows.add(cur);
//                    continue;
//                }
//
//                if (line.contains("START")) {
//                    long start = System.currentTimeMillis();
//                    cur = new MSTwoWindow();
//                    line = br.readLine();
//                    String[] range = line.split("\\s+");
//                    cur.mzLow = Double.parseDouble(range[0]);
//                    cur.mzHigh = Double.parseDouble(range[1]);
//
//                    lineType++;
//                    long time = System.currentTimeMillis() - start;
//                    System.out.println("------------------" + lineType);
//
//                    System.out.println(time);
//
//                    continue;
//                }
//
//            }
//            br.close();
//            freader.close();
//        }
//
//    }

    public int FindWindowWithMZ(double mz, int charge) {
        int len = windows.size();
        int found = FindWindowIndex(mz, len - 1, 0);
        if (found >= 0) {
            return found;
        }
        return -1;
    }

    // given a precursor peptide mass, find the correct isolation window within collection)
    public int FindWindow(double mass, int charge) {
        double mz = MassToMz(mass, charge);
        int len = windows.size();
        int found = FindWindowIndex(mz, len - 1, 0);
        if (found >= 0) {
            return found;
        }
        return -1;
    }
    public static double MassToMz(double mass, int z) {
        return mass/z + AminoAcid.ProtonMass;
    }

    // given mz, return correct isolation window index
    public int FindWindowIndex(double mz, int r, int l) {
        if (r >= l) {
            int mid = l + (r - l) / 2;
            // If the element is present at the
            // middle itself
            if (windows.get(mid).mzHigh >= mz && windows.get(mid).mzLow <= mz) {
                return mid;
            }

            // If element is smaller than mid, then
            // it can only be present in left subarray
            if (windows.get(mid).mzLow > mz) {
                return FindWindowIndex(mz, mid - 1, l);
            }

            // Else the element can only be present
            // in right subarray
            if (windows.get(mid).mzHigh < mz) {
                return FindWindowIndex(mz, r, mid + 1);
            }
        }

        return -1;
    }

//java -Xmx96g -cp  /data/waterlooms-1.0-SNAPSHOT-jar-with-dependencies.jar edu.uw.waterlooms.ReadTrailToFile /data/mouse/ Fig4_mouse_cerebellum_MHRM_R02_T0.mzXML R02Top10Wdecoypsm.csv R02DecoyHtml/
    public static void main(String[] args) throws IOException {
        ReadTrailToFile tst = new ReadTrailToFile();
        tst.OutputMS2TrailInfo(args[0],args[1],args[2],args[3]);
    }
}
