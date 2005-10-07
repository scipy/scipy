""" General purpose classes for plotting utility.

    Many of these classes and function are helpful for laying out the 
    location of element of a graph (titles, tick labels, ticks, etc).
    A simple class for handling properties of graph elements is also
    included.  Everything here is independent of any specific GUI or
    drawing platform
    .
    The classes and functions break out into 3 main sections:
    
    General layout:
        class box_object   -- Helps align rectangular graph elements 
                              relative to one another.  Currently only
                              Int values are allowed for box dimensions.
        class point_object -- Alignment of points with one-another or 
                              box_objects.  
    Property Management:                          
        class poperty_object -- base class for any graph object that 
                                manages a set of properties.  Pretty
                                low-tech compared to graphite and Zope.
                                It could use additions for error checking.
                                
    Axis layout functions:
        auto_ticks          -- find set of tick locations for an axis.
        format_tick_labels  -- convert the numerical tick values to strings
                               for labels
        auto_interval       -- determine an appropriate tick interval for axis
        auto_bounds         -- determine appropriate end points for axis 
        and a few others...
"""


#----------------------------------------------------------
#               General layout classes
#----------------------------------------------------------
from Numeric import *
from scipy_base.fastumath import *
import scipy_base.limits as limits
import scipy_base as misc

LEFT,RIGHT,TOP,BOTTOM = 0,1,2,3 # used by same_as() method

TYPE = Int # all values in a box_object are forced to be integers
           # might change this later when we've looked at round off
           # issues more thoroughly.

class box_object:
    """ Helpful for laying out rectangles.
    
        This class has three general areas of functionality:
            1. Retreive information about a box such as it's size, 
               the x-coordinate of its left side, and others.
            2. Specify the location of one box (rectangle) in relation to 
               another using methods such as left_of(), above(), 
               or center_on().
            3. Change the size (and location) of a box by either choping off
               a portion of the box using trim_right(), trim_left(), and 
               others or inflating the box by a percentage.               
    """
    def __init__(self,topleft,size):
        self.topleft = asarray(topleft,TYPE)
        self._size = asarray(size,TYPE)

    #--------------- interrogation functions --------------------#
    def left(self): 
        "Return the x-coordinate of the left edge."
        return self.topleft[0]
    def right(self): 
        "Return the x-coordinate of the right edge."
        return self.topleft[0] + self.width()
    def top(self): 
        "Return the y-coordinate of the top edge."
        return self.topleft[1]
    def bottom(self):
        "Return the y-coordinate of the bottom edge." 
        return self.topleft[1] + self.height()
    def center_x(self): 
        "Return the x-coordinate of the boxes center."
        return self.topleft[0] + .5 * self.width()
    def center_y(self): 
        "Return the y-coordinate of the left edge."
        return self.topleft[1] + .5 * self.height()
    def size(self):
        "Return the size of the box as array((width,height))." 
        return self._size             
    def width(self):
        "Return the width of the box (length along x-axis)"
        return self.size()[0]    
    def height(self):
        "Return the width of the box (length along y-axis)"
        return self.size()[1]    
    def center(self):
        "Return the coordinates of the boxes center as an array"
        return self.topleft + .5 * self.size()
    def radial(self,angle):
        """ Length of a ray extending from the center of the box at the
            specified polar angle (in radians) to the edge of the box"""
        # Could have some div_by_zero errors if h or w = 0
        #   --> haven't thought about this much
        w,h = self.size()
        slope = abs(tan(angle))
        if slope > float(h)/w:
            s1 = h*.5
            s2 = s1/slope           
        else:            
            s1 = w*.5
            s2 = s1*slope            
        return sqrt(s1**2+s2**2)

    #------------- Moving and resizing methods ----------------#                
    def move(self,location):
        """ Move box so that its top left corner is located at
            the specified location.  Location is an (x,y) pair."""
        self.topleft = asarray(location,TYPE)        
    def translate(self,offset):
        """ Shift the box's location by the amount specifed in
            offset.  Offset is an (x,y) pair."""
        self.topleft = self.topleft + asarray(offset,TYPE)
    def set_size(self,sz):
        """ Set the size of the box. sz is an (x,y) pair."""
        self._size = asarray(sz,TYPE)
    def set_width(self,width):
        """ Set the width of the box. sz is a numeric value"""
        self.set_size((width,self.height()))
    def set_height(self,height):
        """ Set the height of the box. sz is a numeric value."""
        self.set_size((self.width(),height))

    #-------------------- Relational Positioning ----------------#        
    def above(self,other_box,margin=0):
        """Move box so that its bottom edge has same y-coord of the
           top edge of other_box.  The size does not change.  Optionally,
           margin specifies the gap between the boxes."""
        self.topleft[1] = other_box.top() - self.height() - margin
    def below(self,other_box,margin=0):
        """Move box so that its top edge has same y-coord of the 
           bottom edge of other_box.  The size does not change.  Optionally,
           margin specifies the gap between the boxes."""
        self.topleft[1] = other_box.bottom() + margin
    def left_of(self,other_box,margin=0):
        """Move box so that its right edge has same x-coord of the
           left edge of other_box.  The size does not change.  Optionally, 
           margin specifies the gap between the boxes."""
        self.topleft[0] = other_box.left() - self.width() - margin
    def right_of(self,other_box,margin=0):
        """Move box so that its left edge has same x-coord of the 
           right edge of other_box.  The size does not change.  Optionally,
           margin specifies the gap between the boxes."""
        self.topleft[0] = other_box.right() + margin
    def center_on_x_of(self,other_box): 
        """Move box so that x-coord of its center equals the x-coord of 
           other_box's center. The size does not change."""
        self.topleft[0] = other_box.center_x() - .5 * self.width()            
    def center_on_y_of(self,other_box):
        """Move box so that y-coord of its center equals the y-coord of 
           other_box's center. The size does not change."""
        self.topleft[1] = other_box.center_y() - .5 * self.height()            
    def center_on(self,other_box):
        """Move box so that its center is the same as other_box's center.
           The size does not change."""
        self.center_on_x_of(other_box)
        self.center_on_y_of(other_box)
    def same_as(self,other_box,edge,margin=0):
        """Set the specified edge of this box to be aligned with the
           same edge of other_box.  The size does not change.
        """   
        if edge == LEFT: 
            self.topleft[0] = other_box.left() + margin
        elif edge == RIGHT: 
            self.topleft[0] = other_box.right() - margin
        elif edge == TOP: 
            self.topleft[1] = other_box.top() + margin
        elif edge == BOTTOM: 
            self.topleft[1] = other_box.bottom() - margin                    
        else:
            raise ValueError, "edge should only be plt.xxx where xxx is" \
                              " LEFT,RIGHT,TOP, or BOTTOM"                              
    def radial_offset_from(self,other_box,angle,margin=0):
        """ Move this box so that its center is aligned along a radial ray
            out of other_box's center at the specified angle (in radians),
            and its edge touches the appropriate edge of other_box.
            Optionally, margin specifies the size of the gap (along the ray)
            between the two boxes.
        """ 
        dist = other_box.radial(angle) + margin + self.radial(angle)
        new_center = other_box.center() + rotate((dist,0),(0,0),angle)
        self.center_on(point_object(new_center.astype(TYPE)))
                
    #------------------ Trim the edges of the box ----------------#                 
    def trim_top(self,amount,margin=0):
        """ Remove amount from the top edge of the box.  The size and
            topleft corner of the box are altered appropriately.  Margin
            specifies an additional amount to be removed.
        """     
        self.topleft[1] = self.topleft[1] + amount
        self.set_height( self.height() - amount - margin)
    def trim_left(self,amount,margin=0):    
        """ Remove amount from the left edge of the box.  The size and
            topleft corner of the box are altered appropriately.  Margin
            specifies an additional amount to be removed.
        """     
        self.topleft[0] = self.topleft[0] + amount
        self.set_width( self.width() - amount - margin)
    def trim_right(self,amount,margin=0):    
        """ Remove amount from the right edge of the box.  The size and
            is altered appropriately.  Margin specifies an additional amount
            to be removed.
        """     
        self.set_width(self.width() - amount - margin)
    def trim_bottom(self,amount,margin=0):    
        """ Remove amount from the bottom edge of the box.  The size 
            is altered appropriately.  Margin specifies an additional 
            amount to be removed.
        """             
        self.set_height(self.height() - amount - margin)
    def trim_all(self,amount,margin=0):
        """ Remove amount from the all edges of the box.  The size and
            topleft corner of the bow are altered appropriately.  Margin
            specifies an additional amount to be removed.
        """     
        self.trim_top(amount,margin)
        self.trim_bottom(amount,margin)
        self.trim_left(amount,margin)
        self.trim_bottom(amount,margin)
    
    def inflate(self,percent):
        """ Expand the box in all directions by the specified percentage.
            Values < 1 shrink the box.  The center of the box remains the same.
        """
        old_size = self.size()
        inflated = old_size * percent
        shift = .5 * (old_size - inflated)
        self.topleft = floor(self.topleft + shift).astype(TYPE)
        self.set_size(floor(inflated).astype(TYPE))
    
    #------------------ test point relationships ----------------#
    def contains(self,pos):
        l,r,t,b = self.left(),self.right(),self.top(),self.bottom()
        #print self.text, self.size(), l,r,t,b
        #print pos[0],pos[1]
        if ((l <= pos[0] <= r) and (t <= pos[1] <= b)):
            return 1
        return 0
        
class point_object(box_object):
    """ Useful for laying points (and boxes) in relation to one another.
        
        point_object is a degenerate case of a box_object that has
        zero size.  The set_size() method is not allowed for points. As
        a result, size altering methods such as set_height and trimming
        functions are not valid.  Also, the radial() method always returns
        zero.  Other than that, the standard box_object methods work.
        point_objects can also be used as arguments to box_object methods.
    """
    def __init__(self,pt):
        box_object.__init__(self,pt,(0,0))        
    def radial(self,angle):        
        return 0        
    def set_size(self,sz):
        # this'll catch all set_size, trim_xxx and other inappropriate methods
        raise TypeError, "can't set the size of a point_object"

def bounding_points(objects):
    l = min(map(lambda x: x.left(),objects))
    t = min(map(lambda x: x.top(),objects))
    r = min(map(lambda x: x.right(),objects))
    b = min(map(lambda x: x.bottom(),objects))
    return (l,t),(r,b)

def bounding_box(objects):
    l = min(map(lambda x: x.left(),objects))
    t = min(map(lambda x: x.top(),objects))
    r = min(map(lambda x: x.right(),objects))
    b = min(map(lambda x: x.bottom(),objects))
    return box_object((l,t),(r-l,b-t))

#----------------------------------------------------------
#               Simple Property class
#----------------------------------------------------------

class property_object:
    """ Base class for graph object with properties like 'color', 'font', etc.
    
        Many object have a standard set of properties that describe the
        object.  Text objects, for example, have a color, font, and style.
        These attirbutes often have default values and perhaps only a range
        of acceptable values.  Every class derived from this one must have
        the variable "_attributes" defined at the class level.  The keys of
        this dictionary are the names of the class attributes.  The value
        for each key is a list with three entries.  The first is the default
        value for the attribute, the second is a list of acceptable values
        and the third is a short text description of the attribute.  FOr
        example:
          class text_object(property_class):
            _attributes = { 'color': ['black',['black','red','blue'],
                                      "The color of the text" ],
                            'style': ['normal',['normal','bold'],
                                      "The style of the text" ],}

        Currently only the first entry is used.  The second is rather limited
        in functionality, but might be useful for automatically generating
        dialog boxes.
        
        Graphite has a more flexible property system, but it is somewhat
        complex.  Zope also has a nice property structure.  If we run into
        major limitations with this 15 line implementation, we might look
        here for inspiration (or adoption...).
        
        There is no type safety enforced here which could cause so problems
        for unsuspecting users.  Maybe the optional typing of future Python
        versions will mitigate this.
    """    
    def __init__(self, kw, **attr):
        attrs = {};attrs.update(kw);
        if attr: attrs.update(attr)
        for name, value in self._attributes.items():
            val = value[0]
            try:
                val = attrs[name]
            except KeyError: pass
            self.__dict__[name] = val
    def reset_default(self):
        """ Reset the objects attributes to their default values.
        """
        for name, value in self._attributes.items():
            self.__dict__[name] = value[0]

    def clone_properties(self,other):
            """ Reset the objects attributes to their default values.
            """
            for name in other._attributes.keys():
                self.__dict__[name] = other.__dict__[name]

#----------------------------------------------------------#
#-------------- Axis and Tick mark utilities --------------#
#----------------------------------------------------------#

def rotate(pts,center,angle):
    """ pts is a Nx2 array of N (x,y) point pairs.  The points
        are rotated counter-clockwise around a center (x,y) point
        the specified angle (in radians).
    """
    t_pts = asarray(pts) - asarray(center)
    transform = array( ((cos(angle),-sin(angle)),
                        (sin(angle), cos(angle))) )
    rotated = transpose(dot(transform, transpose(t_pts)))
    t_pts = around(rotated + center)
    return t_pts.astype(TYPE)

def log2(num):
    """ Log base 2 of a number ( or array)
        
        !! 1e-16 is here to prevent errors when log is 0
    """
    if is_number(num) and num == 0: 
        num = num + 1e-16
    elif type(num) == type(array([])):
        putmask(num,equal(num,0.),1e-16)    
    return log10(num)/log10(2)

def is_base2(range):
    " True if value is base 2 (2, 4, 8, 16, ...)"
    l = log2(range)
    return (l == floor(l) and l > 0.0)

def is_number(val):
    " Returns 1 if value is an integer or double, 0 otherwise."
    return (type(val) in [type(0.0),type(0)]) 
    
default_bounds = ['auto','auto','auto']
def auto_ticks(data_bounds, bounds_info = default_bounds):
    """ Find locations for axis tick marks.
    
        Calculate the location for tick marks on an axis. data_bounds is a 
        sequence of 2 numbers specifying the maximum and minimum values of 
        the data along this axis. bounds_info is a sequence of 3 values that
        specify how the axis end points and tick interval are calculated. An
        array of tick mark locations is returned from the function.  The
        first and last tick entries are the axis end points.
        
        data_bounds -- (lower,upper).  The maximum and minimum values of the
                       data long this axis.  If any of the settings in 
                       bounds_info are 'auto' or 'fit', the axis properties
                       are calculated automatically from these settings.
        bounds_info -- (lower,upper,interval).  Each entry can either be
                       a numerical value or a string.  If a number,the axis 
                       property is set to that value.  If the entry
                       is 'auto', the property is calculated automatically.
                       lower and upper can also be 'fit' in which case
                       the axis end points are set equal to the values
                       in data_bounds.
    """
    # pretty ugly code...
    # man, this needs some testing.
    
    # hmmm, all the nan stuff should have been handled before this 
    # point...
    #from misc import nan_to_num, isinf
    #if is_number(bounds_info[0]): lower = nan_to_num(bounds_info[0])
    #else:                         lower = nan_to_num(data_bounds[0])
    #if is_number(bounds_info[1]): upper = nan_to_num(bounds_info[1])
    #else:                         upper = nan_to_num(data_bounds[1])
    #if is_number(bounds_info[2]): interval = nan_to_num(bounds_info[2])
    #else:                         interval = bounds_info[2]

    if is_number(bounds_info[0]): lower = bounds_info[0]
    else:                         lower = data_bounds[0]
    if is_number(bounds_info[1]): upper = bounds_info[1]
    else:                         upper = data_bounds[1]
    interval = bounds_info[2]

    #print 'raw input:', lower,upper,interval
    #print 'raw interval:', interval       
    if interval in ['linear','auto']:
        rng = abs(upper - lower)
        if rng == 0.:
            # anything more intelligent to do here?
            interval = .5
            lower,upper = data_bounds + array((-.5,.5))
        if is_base2(rng) and is_base2(upper) and rng > 4:
            if rng == 2: 
                interval = 1
            elif rng == 4: 
                interval = 4
            else: 
                interval = rng / 4 # maybe we want it 8
        else:
            interval = auto_interval((lower,upper))
    elif type(interval) in [type(0.0),type(0)]: 
        pass   
    else:
        #print 'interval: ', interval
        raise ValueError, interval + " is an unknown value for interval: " \
                          "  expects 'auto' or 'linear', or a number"
    
    # If the lower or upper bound are set to 'auto', 
    # calculate them based on the newly chosen interval.
    #print 'interval:', interval
    auto_lower,auto_upper = auto_bounds(data_bounds,interval)
    if bounds_info[0] == 'auto':
        lower = auto_lower
    if bounds_info[1] == 'auto':
        upper = auto_upper        
    
    # again, we shouldn't need this
    # cluge to handle inf values    
    #if isinf(lower):
    #    lower = nan_to_num(lower) / 10
    #if isinf(upper):
    #    upper = nan_to_num(upper) / 10    
    #if isinf(interval):
    #    interval = nan_to_num(interval) / 10    
        # if the lower and upper bound span 0, make sure ticks
    # will hit exactly on zero.
    
    if lower < 0 and upper > 0:
        #print 'arrrgh',0,upper+interval,interval
        hi_ticks = arange(0,upper+interval,interval)
        low_ticks = - arange(interval,-lower+interval,interval)        
        ticks = concatenate((low_ticks[::-1],hi_ticks))
    else:
        # othersize the ticks start and end on the lower and 
        # upper values.
        ticks = arange(lower,upper+interval,interval)
    #cluge
    if len(ticks) < 2:
        ticks = array(((lower-lower*1e-7),lower))
    #print 'ticks:',ticks
    if bounds_info[0] == 'fit': ticks[0] = lower
    if bounds_info[1] == 'fit': ticks[-1] = upper
    return ticks


def format_tick_labels(ticks):
    """ Convert tick values to formatted strings.
        Definitely needs some work
    """
    return map(str,ticks)        
    #print ticks
    #if equal(ticks,ticks.astype('l')):
    #    return map(str,ticks.astype('l'))
    #else:
    #    return map(str,ticks)        
    
        

def calc_bound(end_point,interval,end):
    """ Find an axis end point that includes the the value end_point.  If the
        tick mark interval results in a tick mark hitting directly on the 
        end_point, end_point is returned.  Otherwise, the location of the tick
        mark just past the end_point is returned. end is 'lower' or 'upper' to
        specify whether end_point is at the lower or upper end of the axis.
    """
    quotient,remainder = divmod(end_point,interval)    
    if not remainder: return end_point

    c1 = axis_bound = (quotient + 1) * interval
    c2 = axis_bound = (quotient) * interval
    if end == 'upper': return max(c1,c2)
    if end == 'lower': return min(c1,c2)                

    
def auto_bounds(data_bounds,interval):
    """ Calculate an appropriate upper and lower bounds for the axis from
        the the data_bounds (lower, upper) and the given axis interval.  The
        boundaries will either hit exactly on the lower and upper values
        or on the tick mark just beyond the lower and upper values.
    """    
    data_lower,data_upper = data_bounds
    lower = calc_bound(data_lower,interval,'lower')
    upper = calc_bound(data_upper,interval,'upper')
    return array((lower,upper))

def auto_interval(data_bounds):
    """ Calculate the tick intervals for a graph axis.

        Description:        
        The boundaries for the data to be plotted on the axis are:
            data_bounds = (lower,upper)
                
        A choice is made between 3 to 9 ticks marks (including end points)
        and tick intervals at 1, 2, 2.5, 5, 10, 20, ...
        
        Returns:
        interval -- float. tick mark interval for axis            
    """    
    range = float(data_bounds[1]) - float(data_bounds[0])
    # We'll choose from between 2 and 8 tick marks
    # Favortism is given to more ticks:
    #   Note reverse order and see CLUGE below...
    divisions = arange(8,2.,-1.) #(7,6,...,3)
    # Calculate the intervals for the divisions
    candidate_intervals = range / divisions
    # Get magnitudes and mantissas for each candidate
    magnitudes = 10.**floor(log10(candidate_intervals))
    mantissas = candidate_intervals / magnitudes
    
    # list of "pleasing" intervals between ticks on graph
    # only first magnitude listed, higher mags others inferred
    magic_intervals = array((1.,2.,2.5,5.,10.))
    # calculate the absolute differences between the candidates
    # (with mag removed) and the magic intervals
    differences = abs(magic_intervals[:,NewAxis] - mantissas) 

    # Find the division and magic interval combo
    # that produce the smallest differences.
    # also:
    # CLUGE: argsort doesn't preserve the order of
    # equal values, so we subtract a small , index
    # dependent amount from each difference to
    # force correct ordering
    sh = shape(differences)
    small = 2.2e-16 *arange(sh[1]) * arange(sh[0])[:,NewAxis]
    small = small[::-1,::-1] #reverse order
    differences = differences - small
    # ? Numeric should allow keyword "axis" ? comment out for now
    #best_mantissa = minimum.reduce(differences,axis=0)    
    #best_magic = minimum.reduce(differences,axis=-1)
    best_mantissa = minimum.reduce(differences,0)    
    best_magic = minimum.reduce(differences,-1)
    magic_index = argsort(best_magic)[0]
    mantissa_index = argsort(best_mantissa)[0]
    # The best interval is the magic_interval
    # multiplied by the magnitude of the best mantissa
    interval = magic_intervals[magic_index]
    magnitude = magnitudes[mantissa_index]
    #print differences
    #print 'results:', magic_index, mantissa_index,interval, magnitude
    #print 'returned:',interval*magnitude
    result = interval*magnitude
    if result == 0.0:
        result = limits.float_epsilon
    return result
