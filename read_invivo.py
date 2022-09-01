from urllib.request import urlopen
import sys
import time

fname = 'invivodec02'
delim = '|'

dataset = []

columns = ['tumor', 'strain', 'site', 'schedule', 'route', 'vehicle', 'parameter', 'NSC', 'unit', 'dose', 'T/C', 'endpoint', 'survived', 'x' ]

def read_dict(fname):
	""" read mapping dictionaries from the NCI data instructions """
	result = dict()
	with open(fname, 'r') as file:
		for line in file.readlines():
			key = line[0:2].strip()
			data = line[2:].strip()
			result[key] = data

	return result


""" map the data """
model = read_dict('screening_models.txt')
strain = read_dict('strain.txt')
site  = read_dict('site.txt')
schedule = str
route  =   str
vehicle = read_dict('vehicles.txt')
param = read_dict('parameter.txt')
compound = int
unit =  str
dose =  float
tc = int
endpoint = int
survived = str
x = str
mappers = [model, strain, site, schedule, route, vehicle, param, compound, unit, dose, tc, endpoint, survived, x ]



def mapdata(i, value):

	func = mappers[i] 

	if type(func) is dict:
		return func[value] 
	elif func is None:
		return value
	else:
		try:
			return func(value)
		except:
			return value


def duration(schedule):
	#Q01DX005
	try:
		every = int(schedule[1:3])
	except:
		every = 0

	time  = schedule[3]
	try:
		count  = int(schedule[-3:])
	except:
		count = 0

	duration = float(every * count)
	# convert to days
	if time == 'H' : duration /= 24.0
	if time == 'M' : duration /= 24.0 * 60
	return duration


def assess(datamap):
	readout = datamap['parameter']
	result =   datamap['T/C']
	assess = False
	if 'survival time' in readout.lower():
		if result > 130: assess = True
		datamap['active'] = assess
	if 'tumor weight' in readout.lower():
		if result < 20: assess = True
		datamap['active'] = assess
	
	return datamap


def add_tags(sdfile, datamap):
	sdfile = sdfile[:sdfile.find('$$$$')]
	for tag in datamap:
		sdfile += '>  <' + tag + '>\n'
		sdfile += str(datamap[tag]) + '\n\n'

	sdfile += '$$$$'
	return sdfile


def createurl(number):
	queryformat = 'sdf'
	queryurl = 'https://cactus.nci.nih.gov/chemical/structure/NSC%d/file?format=%s' % (number, queryformat)
	return queryurl


#------------------------------start processing --------------------------

start = time.time()
with open(fname, 'r') as file:
	for line in file.readlines():
		data = line.split(delim)
		datamap = dict() 
		for i, item in enumerate(data):
			item = item.strip()
			tag  = columns[i]
			item = mapdata(i, item)
			datamap[tag] = item	
	
		if 'schedule' in datamap:	
			datamap['duration'] = duration( datamap['schedule'])

		datamap = assess(datamap)

		dataset.append(datamap)

sys.stderr.write('read %7d records in %6.1f\n' % (len(dataset), (time.time()-start)) )

# cache structures to lower internet traffics
cached_structure = dict()
cache_size = 10000
cache_hits = 0
cache_miss = 0
log_interval = 1000
start = time.time()

# read the SDfiles and annotate
for i, datamap in enumerate(dataset):
	cmpd = int(datamap['NSC'])
	if i > 0 and i % log_interval == 0:
		pct = 100 * i/len(dataset)
		elapsed = time.time() - start
		remain = (i/elapsed * (len(dataset) - i))/60
		sys.stderr.write('wrote  %7d sdfiles %7.5f%% cache_hit %6d cache_miss %6d remaining time %6.1f m\n' % (i, pct, cache_hits, cache_miss, remain) )

	if cmpd in cached_structure:
		structure = cached_structure[cmpd]
		cache_hits += 1
	else:
		cache_miss += 1
		url = createurl(cmpd)
		structure = urlopen(url).read().decode('utf-8')
		cached_structure[cmpd] = structure
		if len(cached_structure) > cache_size:
			for item in cached_structure:
				cached_structure.pop(item)
				break

	structure = add_tags(structure, datamap)
	print(structure)

sys.stderr.write('END')

