"""
"""

from plot_utility import *
from wxPython     import wx
from colormap     import colormap_map

#-----------------------------------------------------------------------#
#------- Attribute list/functions for objects in wxPython --------------#
#-----------------------------------------------------------------------#

# Actually the wxWindows color dialog is used to set most
# color properties, so this can be ignored.  The list of supported
# name colors is much longer than this

colors = ['black','blue','green','red','yellow','cyan','brown','orange',
          'purple','grey','white','light grey']
line_widths = [1,2,3,4,5]
tick_length = [0,1,2,3,4,5,6,7,8,9,10]

line_style_map = {'solid':wx.wxSOLID,'dot dash': wx.wxDOT_DASH,
                 'dash':wx.wxSHORT_DASH,'dot': wx.wxDOT,
                 'long dash':wx.wxLONG_DASH,'transparent':wx.wxTRANSPARENT}   
line_styles = line_style_map.keys()      

fill_style_map = {'solid':wx.wxSOLID,'stipple': wx.wxSTIPPLE,
                 'back hatch':wx.wxBDIAGONAL_HATCH,
                 'diag hatch': wx.wxCROSSDIAG_HATCH,
                 'forward hatch':wx.wxFDIAGONAL_HATCH,
                 'cross hatch': wx.wxCROSS_HATCH,
                 'horz hatch': wx.wxHORIZONTAL_HATCH,
                 'vert hatch': wx.wxVERTICAL_HATCH,
                 'transparent':wx.wxTRANSPARENT
                 }                 
fill_styles = fill_style_map.keys()

image_type_map = { 'jpg': wx.wxBITMAP_TYPE_JPEG,
                   'jpeg': wx.wxBITMAP_TYPE_JPEG,
                   'png': wx.wxBITMAP_TYPE_PNG,
                   'pcx': wx.wxBITMAP_TYPE_PCX,
                   'tif': wx.wxBITMAP_TYPE_TIF,
                   'tiff': wx.wxBITMAP_TYPE_TIF, }
    
def get_color(in_color):
    """ Convert a color name or  rgb sequence to a wxColour
    """
    if type(in_color) == type(''):
        color = wx.wxNamedColour(in_color) 
    else: 
        r,g,b = in_color
        ##color = wx.wxColour(r,g,b) # mod by GAP 26092003
        color = wx.wxColour(int(r),int(g),int(b)) 
    return color

def default_font():
    """ The default font for text objects.
        This is a method instead of a global variable so that
        each text object will get its own font object that can
        change indepedently of all others.
    """
    return wx.wxFont(12, wx.wxSWISS, wx.wxNORMAL, 
                     wx.wxNORMAL,wx.false,"Arial")                


#-----------------------------------------------------------------------#
#------------------------- Drawing utilities ---------------------------#
#-----------------------------------------------------------------------#

def draw_point_list(start,stop,pen,dc):
    """ Draw lines between the points listed in start and stop.
    
        start and stop are 2xN arrays of x,y point coordinates.  N lines
        are drawn, one between each of the start,stop pairs using the
        specified pen.  The lines are drawn on the given device context, dc.    
    """    
    dc.SetPen(pen)
    for i in range(len(start)):
        pt1 = start[i]
        pt2 = stop[i]
        dc.DrawLine(pt1[0],pt1[1],pt2[0],pt2[1])
    dc.SetPen(wx.wxNullPen)

#-----------------------------------------------------------------------#
#-------------------------- text objects -------------------------------#
#-----------------------------------------------------------------------#

class text_object(box_object,property_object):
    """ Text objects can be placed anywhere in the plot window.  They have 
        "font", "color", and "rotate" attributes as well as "visible" 
        attribute that specifies whether the text is displayed.  The font
        attribute is a wxWindows font.  "color" can either be a named
        color (string such as 'black') or a 3-element sequence specifiy
        rgb values (0 to 255).  Rotatation is specified in degrees.  
        Currently only 0 and 90 degree rotation are supported.
        
        Issues:
          -- The font attribute seems to be ignored when drawing
             rotated text. (wxPython)
          -- Need to think about origin when drawing rotated text.
             The current approach assumes the box_object does not
             rotate as the text rotates, so the insertion point changes
             from the topleft to bottomleft corner.  This inconsistency
             may be problematic in the future, but it was convienient
             for laying out axis titles.
          -- There is some ugliness with calculating the size of 
             the object because text needs to know the dc it is drawing
             in to know its size.  As a result, you must call set_dc(dc)
             before calling size.  To keep from having a dc reference
             hanging around, I added a clear_dc() method.  Necessary??   
    """    
    _attributes = {'color': ['black',colors,"Color of line"],
                   'font':  [None,[],"Font of text object"],
                   'rotate': [0,[], "Angular rotation of text in degrees"],
                   'visible': ['yes',[], "'yes' or 'no'. Is text visible"],
                   }
    #think about rotation and topleft
    def __init__(self,text=None,topleft=None,**attr):

        property_object.__init__(self,attr)
        box_object.__init__(self, (0,0), (0,0))
        if text is None: self.text = ''
        else: self.text = text
        if not topleft: self.topleft = array((0,0),TYPE)
        else:           self.topleft = array(topleft,TYPE)
        self.translate((0,0))
        # make sure each object has its own font object.
        if not self.font: self.font = default_font()
        if not self.color: self.color = 'black'

    # ugly hack to get size calculations to work
    def set_dc(self,dc):
        """ Set the device context of the text_object for size calculations.
        """
        self.dc = dc
    def clear_dc(self):
        """ Get rid of the device context.
        """
        #try/except cludge added right before SPIE conference
        #find error and fix.
        try:
            del self.dc
        except:
            pass
            
    def size(self):
        """ Calculate the width and length of the text in pixels.
        
            If "visible" == 'no', size() always returns (0,0) regardless
            of the text size.  If the rotation angle is 90 degrees,
            the width of the text is actually its height on the screen
            and visa-versa.
            
            Issues:
               -- haven't really thought about how to handle non-orthogonal
                  text rotations.
        """
        if self.visible in ['yes','on',1] and self.text:
            if not hasattr(self,'dc'):
                raise ValueError, "no device context to calculate text " \
                                  "size. Call set_dc() first."            
            preset = (self.dc.GetFont() == self.font)
            if not preset: self.dc.SetFont(self.font)
            sz= array(self.dc.GetTextExtent(self.text))
            
            # Commented the following line out since calling this
            # seems to turn of plotting of any text            
            #if not preset: self.dc.SetFont(wx.wxNullFont)
            
            # should do something here to calculate real width
            # this only works for 90 degree rotations
            if self.rotate: sz = sz[::-1]
        else:
            sz = array((0,0))
        return sz
                                                    
    def set_size(self):
        raise ValueError, "Can't set size of text objects"
        
    def draw(self,dc):
        """ Draw the text on the screen.
        
            If the text is rotated by 90 degrees, the "insertion point"
            is the bottom-left corner.  Otherwise the stanard top-left
            corner is used.
        """
        if self.visible in ['yes','on',1] and self.text:
            self.set_dc(dc)
            color = get_color(self.color)
            dc.SetTextForeground(color)
            dc.SetFont(self.font)
            # hmmm.  rotated text isn't being drawn with the correct font
            # also need to think about origin for rotated text
            if self.rotate == 90:
                ##dc.DrawRotatedText(self.text,self.left(),self.bottom(),
                ##                   self.rotate)    # mod by GAP 26092003
                dc.DrawRotatedText(self.text,int(self.left()),int(self.bottom()),
                                   int(self.rotate))
            else:
                dc.DrawText(self.text,self.left(),self.top()) 

            dc.SetPen(wx.wxNullPen)
            dc.SetFont(wx.wxNullFont)
            # hmmm. really should clear the dc...
            self.clear_dc()
            

class text_window(wx.wxWindow,text_object):
    """
    Window is currently size=0,0 window for event handling.  Text is
    NOT drawn in the window.  Maybe will do this later, but not sure
    if it is worth handling drawing th background (does GTK support
    transparent windows?).
    """
    def __init__(self,parent,text=None,topleft=None,**attr):
        self.plot_canvas = parent
        wx.wxWindow.__init__(self,parent,-1)
        #not quite right with attr - 1.6 makes this cleaner
        text_object.__init__(self,text,topleft)
        self.SetSize((0,0))

    def format_popup(self,pos):
        menu = wx.wxMenu()
        menu.Append(500, 'Change Text', 'Change Text')
        wx.EVT_MENU(self, 500, self.OnText)
        menu.Append(600, 'Change Font', 'Change Text Font')
        wx.EVT_MENU(self, 600, self.OnFont)
        #print 'fp:',self.text
        menu.UpdateUI()
        self.PopupMenuXY(menu,pos[0],pos[1])
        
    def OnFont(self,event):
        data = wx.wxFontData()
        data.SetColour(get_color(self.color))
        data.SetInitialFont(self.font)
        dlg = wx.wxFontDialog(self, data)
        if dlg.ShowModal() == wx.wxID_OK:
            data = dlg.GetFontData()
            self.font = data.GetChosenFont()
            color = data.GetColour()
            self.color = color.Red(),color.Green(),color.Blue()
            self.plot_canvas.update()
        dlg.Destroy()

    def OnText(self,event):
        #print 'on text:',self.text
        dlg = wx.wxTextEntryDialog(self.GetParent(), 'New Text',
                            defaultValue=self.text)
        if dlg.ShowModal() == wx.wxID_OK:
            self.text = dlg.GetValue()
            self.plot_canvas.update()     
        dlg.Destroy()

        
#-----------------------------------------------------------------------#
#-------------------------- legend_object ------------------------------#
#-----------------------------------------------------------------------#

class legend_object(property_object,box_object):
    """ Draws a legend to the lines plotted on the graph.
    
        Unimplemented.
    """
    _attributes = {
    'position': ['auto',[],"Specify location for legend - auto only for now"],
    'font':       [None,[],"Font of tick labels"],
    'border':     ['yes',['yes','no'],"Place border around legend"],
    'fill_color': ['white',colors,"Background color in legend"],
    'line_length':[30,[],"Length of lines displayed in legend"],
    'visible':    ['yes',['yes','on','no','off'],"Show the legend"],
    }
    def __init__(self, **attr):
        property_object.__init__(self,attr)
        box_object.__init__(self,(0,0),(0,0))
    def layout(self,graph_lines,graph_area,dc):
        
        self.graph_lines = graph_lines # keep a pointer so we can clone
                                       # properties when drawing
        # create duplicates of the lines            
        
        #figure out where to put it
        #create legend
    def draw(self,dc):
        if self.visible in ['yes','on']:
            # beginning of junk that should be moved
            # in to a line object
            self.labels = []
            self.lines = []
            self.markers = []
            
            for line in self.graph_lines:
                #self.labels.append(text_object(line.name))
                self.labels.append(text_object('bug'))
                self.labels[-1].set_dc(dc)
                #create legend lines and markers and steal properties
                #from actual line.  We'll update the points later
    
            # figure out size
            w,h = 0,0
            margin = 0
            line_pts = []
            line_right = point_object((self.line_length+5,0))
            previous = box_object((0,0),(0,0)) #a 0x0 box at the origin
            for i in range(len(self.labels)):
                label = self.labels[i]
                label.below(previous,margin)
                label.right_of(line_right)
                new_line = poly_line([(0,label.center_y()),
                                      (self.line_length,label.center_y())])
                new_line.clone_properties(line.line)
                self.lines.append(new_line)
                new_marker = poly_marker([(self.line_length/2.,
                                            label.center_y())])
                new_marker.clone_properties(line.markers)
                self.markers.append(new_marker)
                previous = label
    
            for items in self.labels: items.clear_dc()
            # end of junk that should be moved
            for item in self.lines:
                item.draw(dc)
            for item in self.markers:
                item.draw(dc)                
            for item in self.labels:
                item.draw(dc)

#-----------------------------------------------------------------------#
#---------------------------- axis_object ------------------------------#
#-----------------------------------------------------------------------#

class axis_object(property_object):
    """ Class for drawing an axis on a graph.
    
        
    """
    _attributes = {
    'labels_visible':  ['yes',['yes','no'], "Turns on/off labels for axis"],
    'label_color': ['black',colors, "Color of tick labels"],
    'label_font':  [None,[],"Font of tick labels"],
    'label_location':  ['minus',['plus','minus'],"Side of axis for labels"],
    'label_offset':  [4,[],"Offset in pixels from tick marks to tick labels"],
    'ticks_visible': ['yes',['yes','no'], "Turns on/off tick marks for axis"],
    'tick_in':    [5,tick_length, "Tick length in pixels into the graph"],
    'tick_out':   [0,tick_length, "Tick length in pixels out ofthe graph"],
    'tick_weight':[1,line_widths, "Width of tick lines"],
    'tick_color': ['black',colors,"Color of tick labels"],
    'tick_style': ['solid',line_styles,"Tick mark and axis line style"],
    'tick_interval': ['auto',[],"Tick interval: 'auto', 'linear','log', " \
                                "or value"],
    'bounds': [['auto','auto'],[],"Axis boundary: 'auto', 'fit', or value"],
    'grid_visible':  ['yes',['yes','no'], "Turns on and off grid lines" \
                                          " for axis"],
    'grid_color': ['light grey',colors,"Color for grid lines"],
    'grid_weight':[1,line_widths, "Width of grid lines"],
    'grid_style': ['dot',line_styles,"Line style for grid lines"],
                   }
    # axis crossing
    # others...
                    
                    
    def __init__(self, rotate = 0, graph_location='above',**attr):
        # attr to allow or disallow labels        
        property_object.__init__(self,attr)
        self.rotate = rotate
        if graph_location == 'above' or graph_location == 'left':
            self.tick_sign = -1
        elif graph_location == 'below' or graph_location == 'right':
            self.tick_sign = 1                        
        else: 
            raise ValueError, "graph_location can be 'left','right','above" \
                              ", or 'below'. You tried '%s'" % graph_location    
        if not self.label_font:
            self.label_font = default_font()
            self.label_font.SetPointSize(10)
        self.bounds = ['auto', 'auto'] # each object needs its own bounds list
        self.omit_first_label = 0
        self.omit_last_label = 0
    def calculate_ticks(self,data_bounds):
        """ data bounds is (lower bound, upper bound)
           axis_settings is a 3-tuple (lower bound, upper bound, interval)
           Each of its settings can be 'auto', or a numerical value.
           In addition, the bounds can be 'fit' and interval can be 
           'linear' or 'log' (well, not yet).
           'auto' -- bound or interval is calculated for you based on the 
                     data bounds.
           'fit'  -- bound is set to be slightly beyond than data_bounds
            value  -- bound is exactly to specified value.
           Typical settings are the following:
            ['auto','auto','auto'
        """    
        bounds_settings = [self.bounds[0], self.bounds[1], self.tick_interval]
        self.ticks = auto_ticks(data_bounds,bounds_settings)
        # This code tries to omit printing of first and last labels when they 
        # are strange values due to 'fit' settings.
        self.omit_first_label = 0
        self.omit_last_label = 0
        if len(self.ticks) > 3:
            delta = self.ticks[2] - self.ticks[1]
            first_delta = self.ticks[1] - self.ticks[0]
            last_delta = self.ticks[-1] - self.ticks[-2]
            first_diff = abs(first_delta - delta)
            last_diff = abs(last_delta - delta)
            if first_diff > .01*delta:
                self.omit_first_label = 1
            if last_diff > .01*delta:
                self.omit_last_label = 1
            
        self.create_labels()
    
    def create_labels(self):
        self.labels = []
        tick_text = format_tick_labels(self.ticks)
        for text in tick_text:
            label = text_object(text,(0,0),font=self.label_font,
                                           color=self.label_color)
            self.labels.append(label)
        if len(self.labels):    
            if self.omit_first_label:
                self.labels[0] = text_object(' ',(0,0),font=self.label_font,
                                               color=self.label_color)
            if self.omit_last_label:
                self.labels[-1] = text_object(' ',(0,0),font=self.label_font,
                                               color=self.label_color)
        
    
    def max_label_size(self,dc,idx):        
        dim = [0]
        for i in self.labels:
            i.set_dc(dc)
            dim.append(i.size()[idx])
            i.clear_dc()
        return max(dim) + self.tick_out + self.label_offset

    def max_label_width(self,dc):        
        return self.max_label_size(dc,0)
    def max_label_height(self,dc):        
        return self.max_label_size(dc,1)
    def range(self):
        return float(self.ticks[-1]-self.ticks[0])
        
    def layout(self,graph_area,dc):
        if self.rotate == 0:           
            length = graph_area.width()
            grid_length = graph_area.height()
        elif self.rotate == 90 or self.rotate == -90: 
            length = graph_area.height()
            grid_length = graph_area.width()
        else: 
            raise ValueError,'rotate must be 0 or 90. It is %d' % self.rotate            

        # translate to screen units for the ticks.
        scale = length / float(self.ticks[-1] - self.ticks[0])
        tick_locations = floor((self.ticks - self.ticks[0]) * scale)
        
        # make array of points with tick screen units as x value, 0 for y value
        tick_points = []
        for i in tick_locations:
            tick_points.append((i,0))            
        tick_points = array(tick_points)

        # create array of end points for the tick lines
        tick_start =  tick_points + array((0,self.tick_in*self.tick_sign))
        tick_stop =  tick_points - array((0,self.tick_out*self.tick_sign))            
        
        # set up grid lines
        grid_start = tick_points 
        grid_stop =  tick_points + array((0,grid_length*self.tick_sign))
        
        # rotate all of this stuff to the correct axis orientation.        
        angle = -self.rotate*pi/180 # neg. since up is y decreasing in canvas.
        zero_point = array((0,0))
        self.tick_points = rotate(tick_points,zero_point,angle)
        self.tick_start = rotate(tick_start,zero_point,angle)
        self.tick_stop = rotate(tick_stop,zero_point,angle)
        self.grid_start = rotate(grid_start,zero_point,angle)
        self.grid_stop = rotate(grid_stop,zero_point,angle)
        #offset the labels from the tick_points
        if self.label_location == 'plus':
            if tan(self.rotate*pi/180) < 1: label_angle = angle - pi/2
            else: label_angle = angle + pi/2
        elif self.label_location == 'minus':
            if tan(self.rotate*pi/180) < 1: label_angle = angle + pi/2
            else: label_angle = angle - pi/2
        for i in range(len(self.labels)):
            self.labels[i].set_dc(dc)
            pt = point_object(self.tick_points[i])        
            self.labels[i].radial_offset_from(pt,label_angle,
                                              margin=self.label_offset)
            self.labels[i].clear_dc()
        self.labels[i].set_dc(dc)            

        # for fast line drawing.
        # not currently used
        self.single_grid_line = []
        self.single_tick_line = []
        ts = map(tuple,self.tick_start)
        tc = map(tuple,self.tick_points)
        te = map(tuple,self.tick_stop)
        N = len(self.tick_start)
        for i in range(N):
            self.single_tick_line.append(tc[i])
            self.single_tick_line.append(te[i])
            self.single_tick_line.append(ts[i])
            self.single_tick_line.append(tc[i])

        
    def move(self,offset):
        self.tick_start = self.tick_start + offset
        self.tick_stop = self.tick_stop + offset
        self.tick_points = self.tick_points + offset
        self.grid_start = self.grid_start + offset
        self.grid_stop = self.grid_stop + offset
        # need to rebuild fast_line drawing stuff
        for i in self.labels:
            i.translate(offset)
    
    def contains(self,pt,dc):
        mins,maxs = [],[]
        for i in self.labels: i.set_dc(dc)
        lo,hi = bounding_points(self.labels)
        for i in self.labels: i.clear_dc()
        mins.append(lo)
        mins.append(minimum.reduce(self.tick_start))
        mins.append(minimum.reduce(self.tick_stop))
        mins.append(minimum.reduce(self.tick_points))
        maxs.append(hi)
        maxs.append(maximum.reduce(self.tick_start))
        maxs.append(maximum.reduce(self.tick_stop))
        maxs.append(maximum.reduce(self.tick_points))
        lo = minimum.reduce(mins)
        hi = maximum.reduce(maxs)
        sz = hi-lo
        return box_object(lo,sz).contains(pt)
          
    def draw_ticks(self,dc):
        if self.ticks_visible in ['yes','on']:
            style = line_style_map[self.tick_style]
            color = get_color(self.tick_color)
            pen = wx.wxPen(color, self.tick_weight, style)        
            draw_point_list(self.tick_start,self.tick_stop,pen,dc)
            #dc.DrawLines(self.single_tick_line)
        # draw axis line here also
        pt1 = self.tick_points[0]    
        pt2 = self.tick_points[-1]
        draw_point_list([pt1],[pt2],pen,dc)
        
    def draw_grid_lines(self,dc):
        if self.grid_visible in ['yes','on']:
            style = line_style_map[self.grid_style]
            color = get_color(self.grid_color)
            pen = wx.wxPen(color,self.grid_weight, style)        
            draw_point_list(self.grid_start,self.grid_stop,pen,dc)

    def draw_labels(self,dc):
        if self.labels_visible in ['yes','on']:
            for label in self.labels:
                label.draw(dc)
                    
    def draw_lines(self,dc):
        self.draw_grid_lines(dc)
        self.draw_ticks(dc)

    def draw(self,dc):
        self.draw_labels(dc) # really only need to do this sometimes...
        self.draw_lines(dc)

class axis_window(wx.wxWindow,axis_object):
    """
    Window is currently size=0,0 window for event handling.  Text is
    NOT drawn in the window.  Maybe will do this later, but not sure
    if it is worth handling drawing th background (does GTK support
    transparent windows?).
    """
    def __init__(self,parent, rotate = 0, graph_location='above',**attr):
        self.plot_canvas = parent
        wx.wxWindow.__init__(self,parent,-1)
        #not quite right with attr - 1.6 makes this cleaner
        axis_object.__init__(self,rotate,graph_location)
        self.SetSize((0,0))

    def format_popup(self,pos):
        menu = wx.wxMenu()
        menu.Append(600, 'Change Font', 'Change Text Font')
        wx.EVT_MENU(self, 600, self.OnFont)
        menu.UpdateUI()
        self.PopupMenuXY(menu,pos[0],pos[1])
        
    def OnFont(self,event):
        data = wx.wxFontData()
        data.SetColour(get_color(self.label_color))
        data.SetInitialFont(self.label_font)
        dlg = wx.wxFontDialog(self, data)
        if dlg.ShowModal() == wx.wxID_OK:
            data = dlg.GetFontData()
            self.label_font = data.GetChosenFont()
            color = data.GetColour()
            self.label_color = color.Red(),color.Green(),color.Blue()
            self.plot_canvas.update()
        dlg.Destroy()

#-----------------------------------------------------------------------#
#-------------------------- border_object ------------------------------#
#-----------------------------------------------------------------------#

class border_object(property_object):
    _attributes = {
    'ticks_visible': ['yes',['yes','no'], "Turns on/off tick marks for axis"],
    'tick_in':    [5,tick_length, "Tick length in pixels into the graph"],
    'tick_out':   [0,tick_length, "Tick length in pixels out ofthe graph"],
    'tick_weight':[1,line_widths, "Width of tick lines"],
    'tick_color': ['black',colors,"Color of tick labels"],
    'tick_style': ['solid',line_styles,"Tick mark line style"],
    'visible':  ['yes',['yes','no'], "Turns on and off border marks for axis"],
    'color': ['black',colors,"Color for border lines"],
    'weight':[1,line_widths, "Width of border lines"],
    'style': ['solid',line_styles,"Line style for border lines"],
    }
    def __init__(self,**attr):
        # attr to allow or disallow labels        
        property_object.__init__(self,attr)
    
    def add_ticks(self,tick_points, origin, tick_in, tick_out):
        offset = array((0,self.tick_out))
        self.tick_start.append( tick_points + offset + origin )
        offset = array((0,self.tick_in))
        self.tick_stop.append( tick_points + offset + origin )

    def layout(self,graph_area,x_axis,y_axis):
        # top
        x_points = x_axis.tick_points - x_axis.tick_points[0]
        y_points = y_axis.tick_points - y_axis.tick_points[0]
        
        #clear tick lines
        self.tick_start = []
        self.tick_stop = []
        self.tick_center = []
        
        #top
        origin = array((graph_area.left(),graph_area.top()))
        offset = array((0,-self.tick_out))
        self.tick_start.append( x_points + offset + origin )
        offset = array((0,self.tick_in))
        self.tick_stop.append( x_points + offset + origin )
        self.tick_center.append(x_points + origin)
        #right
        origin = array((graph_area.right(),graph_area.bottom()))
        offset = array((self.tick_out,0))
        self.tick_start.append( y_points + offset + origin )
        offset = array((-self.tick_in,0))
        self.tick_stop.append( y_points + offset + origin )
        self.tick_center.append(y_points + origin)
        #bottom
        origin = array((graph_area.left(),graph_area.bottom()))
        offset = array((0,self.tick_out))
        self.tick_start.append( x_points + offset + origin )
        offset = array((0,-self.tick_in))
        self.tick_stop.append( x_points + offset + origin )
        self.tick_center.append(x_points + origin)
        #left
        origin = array((graph_area.left(),graph_area.bottom()))
        offset = array((-self.tick_out,0))
        self.tick_start.append( y_points + offset + origin )
        offset = array((self.tick_in,0))
        self.tick_stop.append( y_points + offset + origin )
        self.tick_center.append(y_points + origin)
        
        self.tick_start = concatenate(self.tick_start)
        self.tick_stop = concatenate(self.tick_stop)
        self.tick_center = concatenate(self.tick_center)
        
        self.border_start = []
        self.border_stop = []
        gr = graph_area
        # could do a rect instead, but it always seems
        # off by 1 or so. Do this for now.
        self.border_start.append((gr.left(),gr.top()))
        self.border_stop.append((gr.left(),gr.bottom()))
        self.border_start.append((gr.left(),gr.top()))
        self.border_stop.append((gr.right(),gr.top()))
        self.border_start.append((gr.left(),gr.bottom()))
        self.border_stop.append((gr.right(),gr.bottom()))
        self.border_start.append((gr.right(),gr.top()))
        self.border_stop.append((gr.right(),gr.bottom()))
        self.border_start = array(self.border_start)
        self.border_stop = array(self.border_stop)
                
        # this stuff is only needed for the fast_draw method
        # needs slight amount of work
        self.single_line = []
        ts = map(tuple,self.tick_start)
        tc = map(tuple,self.tick_center)
        te = map(tuple,self.tick_stop)
        N = len(self.tick_start)
        for i in range(N):
            self.single_line.append(tc[i])
            self.single_line.append(te[i])
            self.single_line.append(ts[i])
            self.single_line.append(tc[i])
            
    def draw(self,dc):
        vis = self.visible in ['yes','on']
        tick_vis = self.ticks_visible in ['yes','on']
        
        #if vis and tick_vis:
        if 0:
            # Really should check styles etc. here to make
            # sure the are the same.
            self.draw_fast(dc)
        else:    
            #draw border
            if vis:
                style = line_style_map[self.style]
                color = get_color(self.color)
                pen = wx.wxPen(color, self.weight, style)        
                draw_point_list(self.border_start,self.border_stop,pen,dc)
            #draw ticks
            if tick_vis:
                style = line_style_map[self.tick_style]
                color = get_color(self.tick_color)
                pen = wx.wxPen(color, self.tick_weight, style)        
                draw_point_list(self.tick_start,self.tick_stop,pen,dc)
        
    def draw_fast(self,dc):
        """ This approach uses wxPythons DrawLines to draw
            the entire border in one call instead of drawing
            each tick individually.  It is 5-10 times faster,
            but draws the border and ticks all in the same style.
            (not a big draw back...)
        """
        style = line_style_map[self.style]
        color = get_color(self.color)
        pen = wx.wxPen(color, self.weight, style)        
        dc.SetPen(pen)
        dc.DrawLines(self.single_line)
        dc.SetPen(wx.wxNullPen)

class poly_points(property_object):
    def __init__(self, points, attr=None):
        property_object.__init__(self,attr)
        self.points = array(points)
        self.scaled = array(self.points,copy=1)
        #self.scaled = map(tuple,self.points)

    def bounding_box(self):
        return minimum.reduce(self.points), \
               maximum.reduce(self.points)

    def scale_and_shift(self, scale=1, shift=0):
        self.scaled = scale*self.points+shift
        #self.scaled = map(tuple,self.scaled)
       
class poly_line(poly_points):
    _attributes = {'color': ['black',colors,"Color of line"],
                   'weight': [1,line_widths, "Weight of Line"],
                   'visible':    ['yes',['yes','no'], "Turn on/off line"],
                   'style': ['solid',line_styles,"solid, dash, dot dash, or dot"]
                  }

    def __init__(self, points, **attr):
        poly_points.__init__(self, points, attr)
    def draw(self, dc):
        if self.visible in ['on','yes']:
            color = get_color(self.color)
            style = line_style_map[self.style]
            dc.SetPen(wx.wxPen(color, self.weight,style))
            try:
                dc.DrawLines(self.scaled)           
            except:
                dc.DrawLines(map(tuple,self.scaled))
            dc.SetPen(wx.wxNullPen)
                    
marker_styles = ['circle','square','dot','triangle','down_triangle',\
                 'cross','plus']                    

class poly_marker(poly_points):
    """ adding a layout method that did most of the computation
        and marker function selection upfront might speed up drawing
    """
    _attributes = {
       'outline_color': ['black',colors,"Color of marker outline"],
       'outline_weight': [1,line_widths, "Weight of Line"],
       'fill_color': ['black',colors,"Fill color of marker"],
       'fill_style': ['solid',fill_styles, "pattern used when filling marker"],
       'size':       [1,[],"Size of marker"],
       'symbol':     ['circle',marker_styles,"Shape used for marker"],
       'visible':    ['yes',['yes','no'], "Turn on/off markers"],       
       }
    def __init__(self,points,**attr):
        poly_points.__init__(self,points,attr)

    def draw(self,dc):
        if self.visible in ['on','yes']:
            color = get_color(self.outline_color)
            weight = self.outline_weight
            size = self.size
            fill_color = get_color(self.fill_color)
            fill_style = fill_style_map[self.fill_style]
            symbol = self.symbol
            dc.SetPen(wx.wxPen(color,weight))
            dc.SetBrush(wx.wxBrush(fill_color,fill_style))
            self._drawmarkers(dc, self.scaled, symbol, size)
            dc.SetPen(wx.wxNullPen)
            dc.SetBrush(wx.wxNullBrush)
    def _drawmarkers(self, dc, coords, marker,size=1):
        f = eval('self._' +marker)
        for xc, yc in coords:
            f(dc, xc, yc, size)

    def _circle(self, dc, xc, yc, size=1):
        ##dc.DrawEllipse(xc-3*size,yc-3*size,6*size,6*size)   #mod by GAP 26092003
        dc.DrawEllipse(int(xc-3*size),int(yc-3*size),6*size,6*size)

    def _dot(self, dc, xc, yc, size=1):
        dc.DrawPoint(xc,yc)

    def _square(self, dc, xc, yc, size=1):
        dc.DrawRectangle(xc-3*size,yc-3*size,6*size,6*size)

    def _triangle(self, dc, xc, yc, size=1):
        dc.DrawPolygon([(-0.5*size*6,0.2886751*size*6),
                        (0.5*size*6,0.2886751*size*6),
                        (0.0,-0.577350*size*6)],xc,yc)

    def _down_triangle(self, dc, xc, yc, size=1):
        dc.DrawPolygon([(-0.5*size*6,-0.2886751*size*6),
                        (0.5*size*6,-0.2886751*size*6),
                        (0.0,0.577350*size*6)],xc,yc)

    def _cross(self, dc, xc, yc, size=1):
        dc.DrawLine(xc-3*size,yc-3*size,xc+3*size,yc+3*size)
        dc.DrawLine(xc-3*size,yc+3*size,xc+3*size,yc-3*size)

    def _plus(self, dc, xc, yc, size=1):
        dc.DrawLine(xc-3*size,yc,xc+3*size,yc)
        dc.DrawLine(xc,yc-3*size,xc,yc+3*size)

class line_object(poly_points):
    """ Combines poly_line and poly_marker into a single
        class
    """
    _attributes = {
       'color': ['auto',['auto','custom'],"'auto' or 'custom'."],
       'marker_type':['auto',['auto','custom'],"'auto' or 'custom'."],
       'line_type': ['auto',['auto','custom'],"'auto' or 'custom'."],
       'name': ['','text string that identifies line'],
       }

    def __init__(self,points,**attr):
        poly_points.__init__(self,points,attr)
        self.markers = poly_marker(points)
        self.line = poly_line(points)
    def clip_box(self,box):
        self.clip = box.left(),box.top(),box.width(),box.height()
    def draw(self,dc):                   
        if hasattr(self,'clip'):
            c = self.clip
            ##dc.SetClippingRegion(c[0]-1,c[1]-1,c[2]+2,c[3]+2) # mod by GAP 26092003
            dc.SetClippingRegion(int(c[0]-1),int(c[1]-1),int(c[2]+2),int(c[3]+2))
        self.line.draw(dc)
        if hasattr(self,'clip'): dc.DestroyClippingRegion()
        if hasattr(self,'clip'):
            c = self.clip
            ##dc.SetClippingRegion(c[0]-5,c[1]-5,c[2]+10,c[3]+10)  # mod by GAP 26092003
            dc.SetClippingRegion(int(c[0]-5),int(c[1]-5),int(c[2]+10),int(c[3]+10))
        self.markers.draw(dc)        
        if hasattr(self,'clip'): dc.DestroyClippingRegion()
        
    def scale_and_shift(self, scale=1, shift=0):
        self.markers.scale_and_shift(scale, shift)
        self.line.scale_and_shift(scale, shift)
    def set_color(self,color):
        self.line.color = color
        self.markers.outline_color = color
        self.markers.fill_color = color


        
class image_object(property_object):
    _attributes = {
      'colormap': ['grey',[],"name of colormap or Nx3 color array"],
      'scale':    ['yes',['yes','no'],"Turn scaling on/off"],
      }
    def __init__(self, matrix,x_bounds=None,y_bounds=None,**attr):
        property_object.__init__(self,attr)        
        if not x_bounds: 
            self.x_bounds = array((0,matrix.shape[1]))
        else:
            # works for both 2 element or N element x
            self.x_bounds = array((x_bounds[0],x_bounds[-1]))
        if not y_bounds: 
            self.y_bounds = array((0,matrix.shape[0]))
        else:
            self.y_bounds = array((y_bounds[0],y_bounds[-1]))
        
        self.matrix = matrix
        self.the_image = self.form_image()

        axis_lengths = array((self.x_bounds[1]-self.x_bounds[0],
                              self.y_bounds[1]-self.y_bounds[0]))
        self.image_pixels_per_axis_unit =array((matrix.shape[1], matrix.shape[0]),Float)/axis_lengths
        self.image_origin = array((self.x_bounds[0],self.y_bounds[0]))
        
    def scale_magnitude(self,image,colormap):
        themin = float(minimum.reduce(ravel(image)))
        themax = float(maximum.reduce(ravel(image)))
        scaled = (image - themin) / (themax-themin) * len(colormap) *.99999
        scaled = scaled.astype('b')
        return scaled
        
    def form_image(self):
        # look up colormap if it si identified by a string
        if type(self.colormap) == type(''):
            try:
                colormap = colormap_map[self.colormap]
            except KeyError:
                raise KeyError, 'Invalid colormap name. Choose from %s' \
                                % `colormap_map.keys()`
        else:
            colormap = self.colormap
        # scale image if we're supposed to.    
        if self.scale in ['yes','on']:
            scaled_mag = self.scale_magnitude(self.matrix,colormap)
        else:    
            scaled_mag = self.matrix.astype('b')
        scaled_mag = clip(scaled_mag,0,len(colormap)-1)
        
        if float(maximum.reduce(ravel(colormap))) == 1.:
            cmap = colormap * 255
        else:
            cmap = colormap    
        
        pixels = take( cmap, scaled_mag)
        del scaled_mag
        # need to transpose pixels in memory...
        bitmap = pixels.astype(UnsignedInt8).tostring()
        image = wx.wxEmptyImage(self.matrix.shape[1],self.matrix.shape[0])
        image.SetData(bitmap)
        return image
        
    def bounding_box(self):
        bb = array((self.x_bounds[0],self.y_bounds[0])), \
             array((self.x_bounds[1],self.y_bounds[1]))
        return bb

    def scale_and_shift(self, scale=1, shift=0,upperleft=None,size=None):
        if scale is 1: scale = array((1,1))
        if shift is 0: shift = array((0,0))
        graph_pixels_per_axis_unit = scale
        self.scale = graph_pixels_per_axis_unit/self.image_pixels_per_axis_unit
        self.origin = shift
        self.graph_upperleft = upperleft
        self.graph_size = size

    def draw(self,dc):
        sz = array((self.the_image.GetWidth(),self.the_image.GetHeight()))
        sz = sz* self.scale
        sz = abs(sz.astype(Int))
        scaled_image = self.the_image.Scale(sz[0],sz[1])
        bitmap = scaled_image.ConvertToBitmap()

        dc.DrawBitmap(bitmap, self.origin[0]+1, 
                      self.origin[1]-scaled_image.GetHeight()+1, wx.false)

import UserList
huge = array((1e308,1e308))
tiny = array((-1e308,-1e308))
class graphic_list(UserList.UserList):
    """ probably should have a layout method
    """
    def __init__(self):
        UserList.UserList.__init__(self)
    
    def bounding_box(self):
        p1 =[]; p2 = []
        for o in self.data:
            p1o, p2o = o.bounding_box()
            p1.append(p1o);p2.append(p2o)
        if len(p1):
            return minimum.reduce(p1), maximum.reduce(p2)
        else: # so big or small that they always loose to other data sets
            return huge,tiny

    def axis_bounds(self):
        p1, p2 = self.bounding_box()        
        return array((p1[0],p2[0])), array((p1[1],p2[1]))
    
    def scale_and_shift(self, scale=1, shift=0):
        for o in self.data:
            o.scale_and_shift(scale, shift)

    def draw(self, canvas):
        for o in self.data:
            o.draw(canvas)

auto_styles = [['blue','circle','solid'], ['green','square','solid'],
               ['red','triangle','solid'],['cyan','down_triangle','solid'],
               ['orange','cross','solid'],['purple','plus','solid'] ]
              
# make a black and white auto?
      
class auto_line_list(graphic_list):
    """ Handles drawing of line_objects with auto styles
        only works with line_objects
    """
    def clip_box(self,box):
        for o in self.data:
            o.clip_box(box)
    def draw(self, canvas):
        idx = 0
        for o in self.data:
            #print 'line vis before:',o.line_type, o.line.visible
            #print 'marker vis before:',o.marker_type, o.markers.visible
            if 'auto' in [o.color,o.line_type,o.marker_type]:
                if idx >= len(auto_styles): 
                    idx = 0
                color,marker,line_type = auto_styles[idx]
                #print 'clm: ',color,marker,line_type
                if 'auto' == o.color: o.set_color(color)
                if 'auto' == o.marker_type and marker:
                    o.markers.symbol = marker
                    o.markers.visible = 'yes'
                if 'auto' == o.line_type and line_type:
                    o.line.style = line_type
                    o.line.visible = 'yes'
                idx = idx + 1
            #print 'line vis after:',o.line_type, o.line.visible
            #print 'marker vis after:',o.marker_type, o.markers.visible
            o.draw(canvas)          
