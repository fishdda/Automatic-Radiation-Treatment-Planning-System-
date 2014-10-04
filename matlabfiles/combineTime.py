cases = [5, 25,104,108,113,123,208,30,79];
fwrt = open('ipoptTime.txt','w');
for c in cases:
	fname = 'ipoptTime_case'+str(c);
	fid = open(fname);
	t = fid.read();
	total = float(t);
	total = total/60;
	fwrt.write(str(total));
	fwrt.write("\n");
	fid.close();
fwrt.close();
