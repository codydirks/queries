from sys import stdout,argv, exit  # , version_info
from urllib2 import urlopen, URLError
from cStringIO import StringIO
import gzip

#from astropy import units as u

from . import SlipyError
from .Framework.Command import Parse, CommandError
from .Framework.Options import Options, OptionsError


class MastError(SlipyError):
    """
    Exception specific to the Simbad module.
    """
pass

def URLEncoded(url):
    """
    URLEncoded( string ):
    Return a string with reserved url characters encoded.
    """
    # check argument type
    if type(url) is not str:
        raise SimbadError('URLEncoded function expects type str!')

    # percent-encoded pair dictionary
    encoded = {
        ' ': '%20',
        '%': '%25',
        '#': '%23',
        '(': '%28',
        ')': '%29',
        '|': '%7c'
        #'+': '%2b'
    }

    return ''.join([
        encoded[character] if character in encoded else character
        for character in list(url)])


def MastScript(instrument,criteria,**kwargs):
	script= [
		'https://archive.stsci.edu/',
		URLEncoded(instrument),
		'/search.php?action=Search&outputformat=CSV&coordformat=dec&',
		URLEncoded(criteria)
		]
	return ''.join(script)

class STISDataset:
    def __init__(self, info_string):
        data=info_string.split(',')
        self.dataset=data[0]
        self.target=data[1]
        self.ra=float(data[2])
        self.dec=float(data[3])
        #self.starttime=float(data[4])
        self.date,self.starttime=data[4].split()
        self.exptime=float(data[5])
        self.grating=data[6]
        self.cenwav=float(data[7])
        if len(data)==9:
            self.angsep=float(data[8])
        else:
            self.angsep=float('nan')

    def __repr__(self):
        return self.dataset+'|'+self.target
    def __str__(self):
        return self.dataset+'|'+self.target

class MastQuery:
    def __init__(self, instrument, criteria, default=float, **kwargs):
        if type(instrument) is not str and type(criteria) is not str:
            raise MastError('Mast.MastQuery function expects str types for arguments.')
        try:
			self.options = Options(kwargs,
                {
                    'parse'  : True    , # parse SIMBAD return file
                    'full'   : False   , # return full line of info
                    'dtype'  : default , # convert return data
                    'is_main': False
				})
			self.parse   = self.options('parse')
			self.full    = self.options('full')
			self.dtype   = self.options('dtype')
			self.is_main = self.options('is_main')
			url=MastScript(instrument,criteria)
			#print url
			response=urlopen(url)
			self.data = str( response.read().decode('utf-8')).strip()

        except OptionsError as err:
            print('\n --> OptionsError:', err.msg )
            raise MastError('Mast.MastQuery was not constructed')

        except URLError as error:
            raise MastError('Failed to contact MAST database')

        if 'not found' in self.data or 'error' in self.data:
            raise MastError('could not be resolved by SIMBAD.')

        if self.parse:
            # pre-parse operation common to all criteria
            self.data = self.data.split('data')[-1]

    def __call__(self):
        """
        Retrieve data from Query
        """
        return self.data

def IUESearch(**kwargs):
	dataset='iue'
	try:
		opts=Options(kwargs,
		{
			'target' : '',
			'ra'     : '',
			'dec'    : '',
			'radius' : '3.0',#radius must be in arcmins
			'cam'    : '3',#defaults is short wav camera only
			'mx'     : 100
		})
		target=opts('target')
		ra=opts('ra')
		dec=opts('dec')
		radius=opts('radius')
		cam=opts('cam')
		mx=str(opts('mx'))
	except OptionsError as err:
		print('\n --> OptionsError:')
		raise SimbadError('Simbad.Query was not constructed')

	critstring=''
	critstring+='iue_cam_no='+str(cam)+'&'
	critstring+='max_records='+mx+'&'
	if target != '':
		critstring += 'target='+target
	elif ra != '' and dec != '':
		critstring += 'ra='+ra+'&dec='+dec
	else:
		raise MastError('Need to Provide RA/DEC or target!')
	critstring += '&radius='+radius
	query=MastQuery(dataset, critstring)
	return [x.split(',') for x in query.data.split('\n')[2:]]

#Expected format of dataset is list entry returned by IUESearch
def GetIUEDataset(dataset):
	#Create url from dataset info
	dataset_name=dataset[0]
	pref=dataset_name[0:3]
	num=dataset_name[3:]
	dataset_url='http://archive.stsci.edu/missions/iue/previews/mx/'+ pref.lower()+'/'+num[0:2]+'000/gz/'+pref[0].lower()+pref[2].lower()+num
	if dataset[7] == 'SMALL':
		dataset_url+='s'
	dataset_url+='.gz'
	#stdout.write(dataset_url+'    \r')
	#stdout.flush()

	#Get data from .gz url
	response=urlopen(dataset_url)
	data=[x.strip() for x in (gzip.GzipFile(fileobj=StringIO(response.read()), mode='r')).read().split('\n')]
	#Parse wavelength info from header
	if data[18][0] != 'w':
		#Low dispersion spectra
		"""
		 Quantity                    Units
	     ^^^^^^^^                    ^^^^^
	     Wavelength                  A
	     Flux                        erg cm-2 sec-1 A-1
	     Standard deviation          erg cm-2 sec-1 A-1
	     Background level            IUE Flux Numbers (FN)
	     Net spectrum                IUE Flux Numbers (FN)
	     IUE data quality (nu) flag  numbered code
	 	"""
		header=data[0:18]
		data=data[18:-1]
		#Uneven numbers of whitespaces, remove
		filtered_data=[]
		for line in data:
			entries=filter(None, line.split(' '))
			filtered_data.append(' '.join(entries))
		data=filtered_data
		data_cols=[[float(x.split(' ')[i]) for x in data] for i in (0,1,2)]#,3,4,5)]
		wav=data_cols[0]
		flux=data_cols[1]
		flux_std_dev=data_cols[2]
	else:
		#High Dispersion spectra
		"""
		 Quantity                    Units
	     ^^^^^^^^                    ^^^^^
	     flux*                       erg cm-2 sec-1 A-1
	     estimated noise**           erg cm-2 sec-1 A-1
	     background level            erg cm-2 sec-1 A-1
		 """
		header=data[0:19]
		data=data[19:-1]
		filtered_data=[]
		for line in data:
			entries=filter(None, line.split(' '))
			filtered_data.append(' '.join(entries))
		data=filtered_data
		data_cols=[[float(x.split(' ')[i]) for x in data] for i in (0,1,2)]
		flux=data_cols[0]
		flux_std_dev=data_cols[1]
		#Generate wavelength array from header info
		wave_string=header[-1]
		wave_start=float(wave_string.split('+')[1].split(',')[0])
		wave_delta=float(wave_string.split('*')[0].split('=')[-1].strip())
		k_max=int(wave_string.split(',')[-1].strip())+1
		wav=[wave_delta*i + wave_start for i in range(0,k_max)]

		#Bad values in high dispersion are flagged as -1 in flux, need to remove
		ln=len(flux)
		i=0
		while i < ln:
			if flux[i] == -1:
				del flux[i]
				del flux_std_dev[i]
				del wav[i]
				i -= 1
				ln -= 1
			i+=1

	return wav,flux,flux_std_dev

def STISSearch(**kwargs):
    """
    Inputs: Series of potential search parameters
    Returns a list of STISDataset objects
    """
    dataset='hst'
    try:
        opts=Options(kwargs,
        {
        	'target'  : '',
        	'ra'      : '',
        	'dec'     : '',
        	'radius'  : '3.0',#radius must be in arcmins
        	'config'  : 'STIS/FUV-MAMA',#defaults is short wav camera only
        	'grating' : 'E140H',
            'obs_type': 'S', # S or C for science or calibration, % for both
            'status'  : '%', # Public or Proprietary, % for both
        	'mx'      : 100
        })
        target=opts('target')
        ra=opts('ra')
        dec=opts('dec')
        radius=opts('radius')
        config=opts('config')
        obs_type=opts('obs_type')
        sci_status=opts('status')
        mx=str(opts('mx'))
        grating=opts('grating')
    except OptionsError as err:
        print('\n --> OptionsError:')
        raise MastError('Mast query was not constructed')
    if grating in ('E140H', 'E140M'):
        config='STIS/FUV-MAMA'
    elif grating in ('E230H','E230M'):
        config='STIS/NUV-MAMA'
    critstring=''
    critstring+='selectedColumnsCSV=sci_data_set_name,sci_targname,sci_ra,sci_dec,sci_start_time,sci_actual_duration,sci_spec_1234,sci_central_wavelength,ang_sep&'
    critstring+='sci_instrume=STIS&sci_instrument_config='+config+'&sci_spec_1234='+grating+'&sci_status='+sci_status+'&sci_aec='+obs_type
    # Error checking search parameters
    if (ra != '' and dec == '') or (ra == '' and dec != ''):
        raise MastError('Need to specify both RA/DEC')
    if (target != '' or ra != '' or dec != '') and radius == '':
        raise MastError('Need to specify search radius.')
    #if radius != '' and (target=='' and ra=='' and dec==''):
    #    raise MastError('Need to specify target or RA/DEC when specifying radius.')

    if target != '':
    	critstring += '&target='+target
    elif ra != '' and dec != '':
    	critstring += '&ra='+ra+'&dec='+dec
    if target != '' and ra != '' and dec != '':
        critstring += '&radius='+radius
    query=MastQuery(dataset, critstring)
    if query.data.strip() != 'no rows found':
        return [STISDataset(x) for x in query.data.split('\n')[2:]]
    else:
        return []
