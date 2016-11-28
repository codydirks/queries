#!/usr/bin/env python
# Copyright (c) Geoffrey Lentner 2015. All Rights Reserved.
#Python 2.7 and functionality extension (c) Cody Dirks 2016
# slipy/SLiPy/Simbad.py
"""
usage: Simbad.py @Attribute <identifier> [**kwargs]

This module allows the user to query the SIMBAD astronomical database from
inside Python or shell commands/scripts.

The 'Attribute' points to a function within this module and indicates
what is to be run. Execute 'Simbad.py @Attribute help' for usage details of
a specific function. Currently available attributes are: `Position`,
`Distance`, and `Sptype`.

The identifier names can be anything recognized by SIMBAD (e.g., Regulus,
"alpha leo", "HD 121475", "del cyg", etc ...) if the name is two parts make
sure to use quotation to enclose it.

The **kwargs is the conventional reference to Python keyword arguments.
These should be specific to the 'Attribute' being pointed to.
"""
from sys import argv, exit  # , version_info
from urllib2 import urlopen, URLError
from string import maketrans, digits

from astropy import units as u

from . import SlipyError
from .Framework.Command import Parse, CommandError
from .Framework.Options import Options, OptionsError


class SimbadError(SlipyError):
    """
    Exception specific to the Simbad module.
    """
pass

#Identifier queries intended to return specific astronomical information
#given an identifier
def IDURLEncoded(url):
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
        '|': '%7c',
        '+': '%2b'
    }

    return ''.join([
        encoded[character] if character in encoded else character
        for character in list(url)])

def IDScript(identifier, criteria):
    """
    Script( criteria ):

    URL script for the SIMBAD astronomical database with added
    URLEncoded(criteria).
    """

    script = [
        'http://simbad.u-strasbg.fr/simbad/sim-script?',
        'script=format%20object%20%22', IDURLEncoded(criteria),
        '%22%0a', IDURLEncoded(identifier)
        ]

    return ''.join(script)

class IDQuery:
    """
    IDQuery( identifier, criteria, **kwargs ):

    Class for querying the SIMBAD astronomical database for 'criteria'
    of 'identifier'.

    kwargs = {
        'parse' : True,  # extract relavent data from SIMBAD return file
        'dtype' : float, # output datatype
    }
    """
    def __init__(self, identifier, criteria, default=float, **kwargs):
        """
        Initiate query to SIMBAD database.
        """
        # check argument types
        if type(identifier) is not str or type(criteria) is not str:
            raise SimbadError('Simbad.Query function expects str'
            'types for arguments.')

        try:
            # keyword argument options for Query
            self.options = Options( kwargs,
                {
                    'parse'  : True    , # parse SIMBAD return file
                    'full'   : False   , # return full line of info
                    'dtype'  : default , # convert return data
                    'is_main': False     # called from Main()
                })

            # assignments
            self.parse   = self.options('parse')
            self.full    = self.options('full')
            self.dtype   = self.options('dtype')
            self.is_main = self.options('is_main')

            # query SIMBAD database
            #with urlopen( Script(identifier, criteria) ) as response:
            #    self.data = str( response.read().decode('utf-8') ).strip()
            response = urlopen( IDScript(identifier, criteria) )
            self.data = str( response.read().decode('utf-8')).strip()


        except OptionsError as err:
            print('\n --> OptionsError:', err.msg )
            raise SimbadError('Simbad.Query was not constructed '
                'for `{}`'.format(identifier))

        except URLError as error:
            raise SimbadError('Failed to contact SIMBAD database for'
            ' `{}`'.format(identifier) )

        if 'not found' in self.data or 'error' in self.data:
            raise SimbadError('`{}` could not be resolved by SIMBAD.'
                .format(identifier))

        if self.parse:
            # pre-parse operation common to all criteria
            self.data = self.data.split('data')[-1]

    def __call__(self):
        """
        Retrieve data from Query
        """
        return self.data

def Position( identifier, **kwargs ):
    """
    Position( identifier, **kwargs ):

    Handle to the Query class with criteria='%C00(d;C)'.
    """
    query = IDQuery( identifier, '%COO(d;C)', **kwargs )

    if query.full:
        query.data = query.data.split('\n')[-1]

    elif query.parse:
        # extract relavent data
        query.data = query.data.split()[-1].split('+')
        if len( query.data ) == 1:
            # dec had '-' not '+'
            query.data    = query.data[0].split('-')
            query.data[1] = '-' + query.data[1]
        # return formatted data type
        query.data = [ query.dtype(pos) * u.degree for pos in query.data ]

    if query.is_main:
        if query.full or not query.parse:
            print( query() )
        else:
            print('{0:.2f} {1:.2f}'.format(*query()))

    else: return query.data

def Distance( identifier, **kwargs ):
    """
    Distance( identifier, **kwargs ):

    Handle to the Query class with criteria='%PLX'
    """
    query =  IDQuery( identifier, '%PLX', **kwargs )

    if query.full:
        data = query.data.split('\n')[-1]

    elif query.parse:

        # extract relavent data
        data = query.data.split()

        if data[1] == '~':
            # nothing found!
            raise SimbadError('No distance found for `{}`'.format(identifier))
        try:

            # convert milli-arcseconds to parsecs
            result = u.pc / ( query.dtype(data[1]) / 1000.0 )

        except ValueError as err:
            raise SimbadError('Use a numeric type for Simbad.Distance!')

        if data[2][0] == '[' and data[2][-1] == ']':
            # TODO: understand SIMBAD error format
            # an uncertainty was given by SIMBAD
            # uncertainty = u.pc * 0.001 * query.dtype(data[1]) * query.dtype(data[2][1:-1]) / (
            # 	0.001 * query.dtype(data[1]) )**2
            # uncertainty = result * query.dtype(data[2][1:-1])
            uncertainty = None

        else: uncertainty = None

    else: data = query.data

    if query.is_main:
        if query.full or not query.parse:
            print( data )
        else:
            print( '{0:.2f}'.format( data ) )

    elif not query.parse:
        return data

    else: return result

    # Measurement(result, error=uncertainty, name='Distance',
    #     notes='Retrieved from SIMBAD database for `{}`'.format(identifier))

def SpType(identifier, **kwargs):
    """
    Handle to the Query class with criteria='%SP'.
    """
    query = IDQuery(identifier, '%SP', **kwargs)

    if query.full:
        # return last full line of query
        query.data = query.data.split('\n')[-1]

    elif query.parse:
        # extract relavent data
        query.data = query.data.split()[1]

    if query.is_main:
        print( query() )

    else: return query()

def IDList(identifier, **kwargs):
    """
    Handle to the Query class with criteria='%IDLIST'.
    With `parse` = True, return a list of alternate IDs for
    the `identifier` provided.
    """
    query = IDQuery(identifier, '%IDLIST', **kwargs)

    if query.parse:
        # extract relavent data
        query.data = query.data.split(':')[-1].strip().split('\n')

    if query.is_main:
        for line in query.data:
            print(line)

    else: return query()

def BVFluxes(identifier, **kwargs):
    query= IDQuery(identifier, '%FLUXLIST(B,V;F,)',**kwargs)
    output=query.data.split('\n')[-1]#.split(',')
    if output == 'simbatch done':
        return [float('NaN'),float('NaN')]
    if query.parse:
        b,v=query.data.split('\n')[-1].split(',')[0:2]
        if b is '':
            bf=float('NaN')
        else:
            bf=float(b)
        if v is '':
            vf=float('NaN')
        else:
            vf=float(v)
        query.data = [bf,vf]
    if query.is_main:
        print query()
    else:
        return query()

def ObjType(identifier, **kwargs):
    query=IDQuery(identifier, '%OTYPE(3)', **kwargs)
    if query.parse:
        query.data=query.data.split('\n')[-1]
    if query.is_main:
        print query()
    else: return query()



#List based queries built around Simbad's criteria searches
#CoordSearch and CritSearch returns lists of objects
def CritURLEncoded(url):
    # check argument type
    if type(url) is not str:
        raise SimbadError('URLEncoded function expects type str!')

    # percent-encoded pair dictionary
    encoded = {
        ' ': '+',
        '%': '%25',
        '#': '%23',
        '(': '%28',
        ')': '%29',
        '|': '%7c',
        '+': '%2b',
        ',': '%2c',
        '&': '%26'
    }
    return ''.join([
        encoded[character] if character in encoded else character
        for character in list(url)])


# Would like to add options to this to choose returned output.
# Being able to select proper motions, parallaxes, etc would be nice
def CritScript(critstring, outputmode='list', mx=100, get_fluxes=True, get_pms=False, get_plx=False):
    mx = int(mx)
    if outputmode.upper() not in ('LIST','COUNT'):
        raise SimbadError('Output mode must be LIST or COUNT!')

    script = [
        'http://simbad.u-strasbg.fr/simbad/sim-sam?',
        'output.format=ASCII&list.idopt=CATLIST&list.idcat=HD', #select HD cat names
        '&list.bibsel=off&list.notesel=off&obj.bibsel=off&obj.notesel=off', #hide bib and notes
        '&coodisp1=[d][2]',#coordinate output format
        '&Criteria=',CritURLEncoded(critstring),
        '&OutputMode=',outputmode,'&maxObject=',str(mx)]
    if get_fluxes:
        script.append('&list.fluxsel=on&U=off&R=off&B=on&V=on') #format list flux/mag display
    else:
        script.append('&list.fluxsel=off')
    if get_pms:
        script.append('&list.pmsel=on')
    else:
        script.append('&list.pmsel=off')
    if get_plx:
        script.append('&list.plxsel=on')
    else:
        script.append('&list.plxsel=off')

    return ''.join(script)

class CritQuery:
        """
        CritQuery( critstring, **kwargs ):

        Class for querying the SIMBAD astronomical database for 'criteria'

        kwargs = {
            'parse' : True,  # extract relavent data from SIMBAD return file
            'dtype' : float, # output datatype
        }
        """
        def __init__(self, criteria, default=float, **kwargs):
            """
            Initiate query to SIMBAD database.
            """
            # check argument types
            if type(criteria) is not str:
                raise SimbadError('Simbad.Query function expects str'
                'types for arguments.')

            try:
                # keyword argument options for Query
                self.options = Options( kwargs,
                    {
                        'parse'  : True    , # parse SIMBAD return file
                        'full'   : False   , # return full line of info
                        'dtype'  : default , # convert return data
                        'is_main': False   , # called from Main()
                        'mode'   : 'LIST'  , # Output mode
                        'mx'     : 100,       # Max number to return
                        'get_coords'   : True,
                        'get_fluxes'   : True,
                        'get_pms'      : False,
                        'get_plx'      : False,
                        'get_spec_type': True
                    })

                # assignments
                self.parse   = self.options('parse')
                self.full    = self.options('full')
                self.dtype   = self.options('dtype')
                self.is_main = self.options('is_main')
                self.mode    = self.options('mode')
                self.mx      = self.options('mx')
                flx          = self.options('get_fluxes')
                pms          = self.options('get_pms')
                plx          = self.options('get_plx')
                # query SIMBAD database
                #with urlopen( Script(identifier, criteria) ) as response:
                #    self.data = str( response.read().decode('utf-8') ).strip()
                url=CritScript(criteria, self.mode, self.mx, flx,pms,plx)
                response = urlopen( url )
                self.data = str( response.read().decode('utf-8')).strip()

            except OptionsError as err:
                print('\n --> OptionsError:', err.msg )
                raise SimbadError('Simbad.Query was not constructed')

            except URLError as error:
                raise SimbadError('Failed to contact SIMBAD database')

            if 'not found' in self.data or 'error' in self.data:
                raise SimbadError('could not be resolved by SIMBAD.')

            if self.parse:
                # pre-parse operation common to all criteria
                self.data = self.data.split('data')[-1]

        def __call__(self):
            """
            Retrieve data from Query
            """
            return self.data

def CoordSearch(lng,lat,rad,**kwargs):
    try:
        opts=Options(kwargs,
            {
                'frame'   : 'icrs',
                'radunit' : 'm',
                'fulldata': False
            })
        frame=opts('frame')
        radunit=opts('radunit')
        fulldata=opts('fulldata')
    except OptionsError as err:
        print('\n --> OptionsError:', err.msg )
        raise SimbadError('Simbad.Query was not constructed')
    except URLError as error:
        raise SimbadError('Failed to contact SIMBAD database for')
    if radunit not in ('d','m','s'):
        raise SimbadError('Unit of radius must be one of d,m,s!')

    #Lat must be preceded by sign
    if lat >= 0:
        lat='+'+str(lat)

    critstring=['region(circle,', frame, ',',
        str(lng),' ',str(lat),',',str(rad),radunit,')']
    query=CritQuery(''.join(critstring), **kwargs)
    if fulldata:
        return query()
    else:
        query.data=query.data.split('\n')
        lst=[]
        for entry in query.data:
            if (len(entry) >0) and (entry[0].isdigit()):
                lst.append(entry.split('|')[1].split('  ')[0])
        query.data=lst
        return query()



class SimbadObject(object):
    def __init__(self,info_string, header):
        data=[x.strip() for x in info_string.split('|')]
        hdrs=[x.strip() for x in header.split('|')]
        coord_string=data[3]
        self.identifier=data[1].split('  ')[0]
        self.objecttype=data[2]
        self.ra=float(coord_string.split(' ')[0])
        self.dec=float(coord_string.split(' ')[1])
        if 'pm' in hdrs and ('~' not in data[hdrs.index('pm')]):
            pms=[x.strip() for x in data[hdrs.index('pm')].split(' ')]
            self.pm_ra = float(pms[0])
            self.pm_dec= float(pms[1])
        else:
            self.pm_ra = float('nan')
            self.pm_dec= float('nan')

        if ('plx' in hdrs) and data[hdrs.index('plx')] != '~':
            self.plx=float(data[hdrs.index('plx')])
        else:
            self.plx=float('nan')

        if 'Mag V' in hdrs:
            mag_start_idx=hdrs.index('Mag U')
            if data[mag_start_idx] == '~':
                self.umag=float('nan')
            else:
                self.umag=float(data[mag_start_idx])
            if data[mag_start_idx+1] == '~':
                self.bmag=float('nan')
            else:
                self.bmag=float(data[mag_start_idx+1])
            if data[mag_start_idx+2] == '~':
                self.vmag=float('nan')
            else:
                self.vmag=float(data[mag_start_idx+2])
            if data[mag_start_idx+3] == '~':
                self.imag=float('nan')
            else:
                self.imag=float(data[mag_start_idx+3])
            if data[mag_start_idx+4] == '~':
                self.rmag=float('nan')
            else:
                self.rmag=float(data[mag_start_idx+4])
        else:
            self.umag=float('nan')
            self.bmag=float('nan')
            self.vmag=float('nan')
            self.imag=float('nan')
            self.rmag=float('nan')

        self.spectraltype=data[-1]

    def __repr__(self):
        return self.identifier
    def __str__(self):
        return self.identifier

#If mode='COUNT', return integer number of hits
#If mode='LIST', returns list of identifiers
def CritSearch(critstring, **kwargs):
    query=CritQuery(critstring, **kwargs)
    try:
        opts=Options(kwargs,
            {
                'mx'       : 100,
                'mode'     :'LIST',
                'full'     : False,
                'get_coords'   : True,
                'get_fluxes'   : True,
                'get_pms'      : False,
                'get_plx'      : False,
                'get_spec_type': True
            })
        mode=opts('mode').upper()
        fulldata=opts('full')
    except URLError as error:
        raise SimbadError('Failed to contact SIMBAD database for')
    if mode=='COUNT':
        return query.data.split('=')[1].strip()
    if fulldata:
        return query()
    else:
        query.data=query.data.split('\n')
        if (len(query.data) == 1):#No results returned
            return []
        if (len(query.data) > 1):# >0 returned
            if query.data[5].split(' ')[0] == 'Number':# >1 returned as list, strip object names from list entries
                header=query.data[7]
                lst=[]
                for entry in query.data:
                    if (len(entry) >0) and (entry[0].isdigit()):
                        #lst.append(entry.split('|')[1].split('  ')[0])
                        lst.append(SimbadObject(entry, header))
                query.data=lst
                return query()
            else:# 1 item return as object(instead of list), need to strip object name from data
                dat=query.data[5].split(' ')
                if dat[1] != 'HD':
                    return [' '.join(dat[1:3])]
                al=maketrans('','')
                nodigs=al.translate(al,digits)
                hdnum=' '.join(dat[1:3]).translate(al,nodigs)
                query.data=['HD '+hdnum]
                return query()

def Main( clargs ):
    """
    Main function. See __doc__ for details.
    """

    if len(clargs) < 2:
        # show usage
        print( __doc__ )
        return 0

    # Available functions for execution
    executable = {
            'Distance' : Distance, # search for parsecs
            'Position' : Position, # search for ra, dec
            'Sptype'   : Sptype  , # search for spectral types
            'IDList'   : IDList    # search for IDs
        }


    try:

        # Parse command line arguments
        function, args, kwargs = Parse( clargs[1:] )

        if not args and not kwargs or args[0] == 'help':
            # show function usage
            print( executable[function].__doc__ )
            return 0

        # run execution
        for identifier in args:
            executable[function]( identifier, is_main=True, **kwargs )

        return 0

    except CommandError as err:
        print(' --> CommandError:', err.msg)
        return 1

    except KeyError as key:
        print(' --> {} was not a recognized function.'.format(key))
        return 1

    except SimbadError as err:
        # don't let uncaught self exception pass if from main.
        print(' --> SimbadError:', err.msg)
        return 1

    except Exception as err:
        print(' --> Unrecognized error with query for `{}`'
                .format(args[0]))
        print(' --> Exception: `{}`'.format(err))
        return 1

if __name__ == '__main__':
    # Main return 0 or 1
    exit( Main( argv ) )
