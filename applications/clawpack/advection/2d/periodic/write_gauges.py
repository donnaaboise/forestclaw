import fclaw_analysis

dim = 2;
gaugedata = fclaw_analysis.GaugeData(dim,min_time_increment=0)

# Format : [id, x, y, t0, t1]
#   id   : integer, identifying the gauge
#  x,y   : Location of the gauge
#  t0,t1 : (t0,t1) interval over which to monitor the gauge.
gaugedata.gauges.append([  1, 0.51, 0.467,   0., 1.e10])

gaugedata.write(data_source='write_gauges.py')
