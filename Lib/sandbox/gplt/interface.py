import new_plot
pyPlot = new_plot
from Numeric import *
from scipy_base.scimath import *


def _intialize():
	global _figure; global _active
	_figure = [pyPlot.Plot()]
	_active = _figure[-1]

_intialize()

def _validate_active():
	global _figure; global _active
	if _active.valid(): return
	try: _figure.remove(_active)
	except ValueError: pass
	try: 
		_active = _figure[-1]		
		_validate_active()
	except IndexError: _intialize()
	
def figure(which_one = None):
	global _figure; global _active
	if which_one is None:
		_figure.append(pyPlot.Plot())
		_active = _figure[-1]
	else:
		try:	_figure.index(which_one)
		except ValueError: _figure.append(which_one)
		_active = which_one

def current():
	_validate_active()
	return _active
	
def close(which_one = None):
	global _figure; global _active
	if which_one is None:
		_active.close()
		try: _active = _figure[-1]
		except IndexError: _intialize()
	elif which_one == 'all':
		for fig in _figure: fig.close()
		_intialize()
	else:
		try: 
			_figure.remove(which_one)
			which_one.close()
		except ValueError:
			which_one.close() 
	if _active == which_one: 
		try: _active = _figure[-1]
		except IndexError: _intialize()

def cmp_imag(x,y):
	if x.real == y.real and x.imag == y.imag: return 0
	if x.real <= y.real: return -1
	if x.real>= y.real: return 1

def zplane(data,title=None, style = 19):
	data = asarray(data,Complex)
	xmin = floor(min(data.real))
	xmax = ceil(max(data.real))
	ymin = floor(min(data.imag))
	ymax = ceil(max(data.imag))
	if(ymax - ymin) == 0:
		ymax = xmax
		ymin = xmin
	if(xmax - xmin) == 0:
		xmax = ymax
		xmin = ymin
	if xmax  < 1.: xmax = 1.
	if ymax  < 1.: ymax = 1.
	if xmin  > -1.: xmin = -1.	
	if ymin  > -1.: ymin = -1.	
	t = arange(100)*2*pi/100.-pi
	_active.plot(sin(t),cos(t),'notitle w l lt 0')
	_active.hold('on')
	_active.plot([xmin,xmax],[0.,0.],'notitle w l lt 0')
	_active.plot([0.,0.],[ymin,ymax],'notitle w l lt 0')
	if title:
		_active.plot(data.real,data.imag, 'title "%s"w p pt %d ps 3' % (title,style))
	else:	
		_active.plot(data.real,data.imag, 'notitle w p pt %d ps 3' % (title,style))
	rr = array(data).tolist()
	rr.sort(cmp_imag)
	SMALL = 1e-2
	for i in rr:
		count = 0
		for j in rr:
			if abs(i-j) < SMALL: count = count + 1
		#print i.real, i.imag, count	
		if count > 1:
			xpos = 1.03 *(( i.real - xmin) / (xmax - xmin))
			ypos = 1.03 *((i.imag - ymin) / (ymax - ymin))
			_active.text((xpos,ypos),"%d"%count)	
	_active.xaxis([xmin,xmax])
	_active.yaxis([ymin,ymax])
	_active._send('set size ratio %f' % ((ymax-ymin)/(xmax-xmin)))	
	_active.hold('off')
	_active.grid('off')
def plot(*data):
	_validate_active()
	apply(_active.plot,data)
def polar(*data):
	_validate_active()
	apply(_active.polar,data)
def autoscale():
	_validate_active()
	_active.autoscale()
def logx(st='on'):
	_validate_active()
	_active.logx(st)
def logy(st='on'):
	_validate_active()
	_active.logy(st)	
def xaxis(rng):
	_validate_active()
	_active.xaxis(rng)
def yaxis(rng):
	_validate_active()
	_active.yaxis(rng)
def zaxis(rng):
	_validate_active()
	_active.zaxis(rng)
def raxis(rng):
	_validate_active()
	_active.raxis(rng)
def taxis(rng):
	_validate_active()
	_active.taxis(rng)
def grid(st):
	_validate_active()
	_active.grid(st)
def hold(st):
	_validate_active()
	_active.hold(st)
def title(t):
	_validate_active()
	_active.title(t)
def xtitle(t):
	_validate_active()
	_active.xtitle(t)
def ytitle(t):
	_validate_active()
	_active.ytitle(t)
def ztitle(t):
	_validate_active()
	_active.ztitle(t)	
def text(xyz,*t):
	_validate_active()
	apply(_active.text,(xyz,)+t )
def label_coord(c):
	_validate_active()
	_active.label_coord(c)
def legend(cmd=''):
	_validate_active()
	_active.legend(cmd)
#3d plotting stuff		
def surf(*data):
	_validate_active()
	apply(_active.surf,data)
def plot3d(*data):
	_validate_active()
	apply(_active.plot3d,data)
def mapping(val):
	_validate_active()
	_active.mapping(val)
def angles(val):
	_validate_active()
	_active.angles(val)
def view(val):
	_validate_active()
	_active.view(val)
def ticlevel(val):
	_validate_active()
	_active.ticlevel(val)
def hidden(val):
	_validate_active()
	_active.hidden(val)
def output(filename,type,options = ''):
	_validate_active()
	_active.output(filename,type,options)

def png(filename):
	_validate_active()
	_active.png(filename)
def save(filename):
	_validate_active()
	_active.png(filename)

def size(size=(1.,1.)):
	_validate_active()
	_active.size(size)

