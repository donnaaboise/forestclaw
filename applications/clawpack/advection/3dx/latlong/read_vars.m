function parms = read_vars()

data = load('setprob.data');

parms.example = data(1);
parms.revs_per_second = data(4);
parms.longitude = data(5:6);
parms.latitude = data(7:8);
parms.maxelev = data(9);
parms.scale = data(10:12);
parms.center = data(13:14);

end