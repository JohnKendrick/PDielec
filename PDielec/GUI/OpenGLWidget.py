import os
import math
import numpy as np
# The following lines seem to fix a problem when running on low end machines
import OpenGL
OpenGL.USE_ACCELERATE = False
OpenGL.ERROR_CHECKING = False
# end of low-end machine fix
from collections     import deque
from OpenGL.GL       import *
from OpenGL.GLU      import *
from OpenGL.GLUT     import *
from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore    import Qt
from PyQt5.QtCore    import QTimer
from PyQt5.QtGui     import QSurfaceFormat
from PDielec.Constants import PI
from PDielec.Utilities import Debug

class OpenGLWidget(QOpenGLWidget):

    format = QSurfaceFormat()
    #format.setDepthBufferSize(24)
    #format.setStencilBufferSize(8)
    format.setSamples(6)
    format.setSwapBehavior(QSurfaceFormat.DoubleBuffer)
    QSurfaceFormat.setDefaultFormat(format)

    def __init__(self, parent, debug=False):
        QOpenGLWidget.__init__(self, parent)
        global debugger
        debugger = Debug(debug,'OpenGLWidget')
        self.debug = debug
        self.viewerTab = parent
        self.notebook = parent.notebook
        self.setMinimumSize(640, 672)
        self.lightingOn = True
        self.diffuseMaterialFactor  = 0.9
        self.ambientMaterialFactor  = 0.4
        self.specularMaterialFactor = 0.7
        self.glintMaterialFactor    = 100.0
        self.diffuseLightFactor     = 0.8
        self.ambientLightFactor     = 0.6
        self.specularLightFactor    = 1.0
        self.linearAttenuation      = True
        self.image_size             = 15.0
        self.white                  = np.array( [1.0, 1.0, 1.0, 1.0] )
        self.spheres                = deque()
        self.cylinders              = deque()
        self.arrows                 = deque()
        self.current_phase          = 0
        self.phase_direction        = 1
        self.number_of_phases       = 1
        self.sphere_slices          = 20
        self.sphere_stacks          = 20
        self.cylinder_slices        = 8
        self.cylinder_stacks        = 2
        self.timer                  = None
        self.show_arrows            = True
        self.timer_interval         = 60
        self.my_width               = None
        self.my_height              = None
        self.background_colour      = None

        self.show_full_screen       = False
        self.writer                 = None
        self.rotation_centre         = np.array( [0.0, 0.0, 0.0] )
        self.matrix =  np.eye( 4, dtype=np.float32)
        self.light_switches = None
        d = 200.0
        self.light_positions = [ [-d, d, d], [ d,-d, d], [-d,-d, d], [ d, d, d], [-d, d,-d], [ d,-d,-d], [-d,-d,-d], [ d, d,-d] ]
        self.lights = [ GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3, GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7 ]

    def enterEvent(self, event):
        debugger.print('enter event')
        self.setFocus()
        self.viewerTab.enterEvent(event)

    def showArrows(self, show):
        debugger.print('show arrows',show)
        self.show_arrows = show

    def timeoutHandler(self):
        self.current_phase += self.phase_direction
        if self.current_phase >= self.number_of_phases:
            self.phase_direction = -1
            self.current_phase = self.number_of_phases -2
        elif self.current_phase < 0:
            self.current_phase = 1
            self.phase_direction = +1
        debugger.print('Timeout - phase', self.current_phase)
        self.update()

    def myMakeCurrent(self):
        self.makeCurrent()
        err = glGetError()
        while err != GL_NO_ERROR:
            err = glGetError()

    def moleculeRotate(self,scale,x,y,z):
        # Rotate molecular frame using glortho
        # Create a 4x4 matrix for the model view
        debugger.print('molecule rotation',x,y,z)
        self.myMakeCurrent()
        glMatrixMode(GL_MODELVIEW)
        glGetFloatv(GL_MODELVIEW_MATRIX, self.matrix)
        glLoadIdentity()
        glRotatef(scale*x, 1.0, 0.0, 0.0)
        glRotatef(scale*y, 0.0, 1.0, 0.0)
        glMultMatrixf(self.matrix)
        self.update()

    def keyPressEvent(self, event):
        self.myMakeCurrent()
        key = event.key()
        modifiers = event.modifiers()
        debugger.print('kepressevent',key,modifiers)
        control = False
        shift = False
        if modifiers & Qt.ShiftModifier:
            shift = True
        if modifiers & Qt.ControlModifier:
            control = True
        if modifiers & Qt.ShiftModifier or modifiers & Qt.ControlModifier:
            amount = 45.0
        else:
            amount = 5.0
        if key == Qt.Key_Left:
            self.moleculeRotate(-amount,0.0,1.0,0.0)
        elif key == Qt.Key_Right:
            self.moleculeRotate(+amount,0.0,1.0,0.0)
        elif key == Qt.Key_Up:
            self.moleculeRotate(-amount,1.0,0.0,0.0)
        elif key == Qt.Key_Down:
            self.moleculeRotate(+amount,1.0,0.0,0.0)
        elif key == Qt.Key_Plus:
            self.zoom(+1.0)
        elif key == Qt.Key_Minus:
            self.zoom(-1.0)
        elif key == Qt.Key_P:
            self.save_movie()
        elif key == Qt.Key_F:
            self.show_full_screen = not self.show_full_screen
            if self.show_full_screen:
                self.showFullScreen()
            else:
                self.showNormal()
        elif key == Qt.Key_Home:
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            self.matrix =  np.eye( 4, dtype=np.float32)
            self.current_phase = int(self.number_of_phases / 2)
            debugger.print('Home key', self.current_phase)
            self.update()

    def save_movie(self, filename):
        import imageio
        # imageio.plugins.ffmpeg.download()
        debugger.print('save_movie', filename)
        if self.timer is not None:
            self.timer.stop()
        writer  = imageio.get_writer(filename, mode='I', fps=24)
        tmpdir  = os.path.dirname(filename)
        tmpfile = os.path.join(tmpdir,'.snapshot.png')
        debugger.print('save_movie', filename)
        for i in range(0,2*self.number_of_phases):
            self.timeoutHandler()
            image = self.grabFramebuffer()
            x = image.width()
            y = image.height()
            modx = x%16
            mody = y%16
            startx = int(modx / 2)
            starty = int(mody / 2)
            x = x - modx + startx
            y = y - mody + starty
            # The image is cropped so that the height and widths are multiples of 16
            image = image.copy(startx, starty, x,y)
            image.save(tmpfile)
            image = imageio.imread(tmpfile)
            writer.append_data(image);
        writer.close()
        os.remove(tmpfile)
        if self.timer is not None:
            self.timer.start()

    def snapshot(self,filename):
        debugger.print('snapshot', filename)
        image = self.grabFramebuffer()
        image.save(filename)

    def translate(self,x,y):
        debugger.print('translate ',x,y)
        self.myMakeCurrent()
        glTranslatef(x, y, 0.0)

    def wheelEvent(self, event):
        debugger.print('Wheel event ' )
        self.myMakeCurrent()
        zoom = event.angleDelta().y()
        self.zoom(zoom)
        self.update()

    def mousePressEvent(self, event):
        self.xAtPress = event.x()
        self.yAtPress = event.y()
        #print('mouse press',self.xAtPress, self.yAtPress)

    def mouseReleaseEvent(self, event):
        self.xAtRelease = event.x()
        self.yAtRelease = event.y()
        #print('mouse release',self.xAtRelease, self.yAtRelease)

    def zoom(self, zoom):
        debugger.print('zoom ', zoom)
        self.myMakeCurrent()
        if zoom > 0:
            zoom_factor = 1.06
        else:
            zoom_factor = 0.94
        debugger.print('zoom factor', zoom_factor)
        glScalef(zoom_factor, zoom_factor, zoom_factor)
        self.update()

    def mouseMoveEvent(self, event):
        self.xAtMove = event.x()
        self.yAtMove = event.y()
        buttons = event.buttons()
        modifiers = event.modifiers()
        if buttons & Qt.LeftButton:
            debugger.print('Mouse event - left button')
            if modifiers & Qt.ShiftModifier or modifiers & Qt.ControlModifier:
                # handle zoom
                xzoom = +1.0*(self.xAtMove - self.xAtPress)
                yzoom = -1.0*(self.yAtMove - self.yAtPress)
                self.zoom(xzoom+yzoom)
            else:
                # handle ordinary rotation
                xrotate = (self.xAtMove - self.xAtPress)*1
                yrotate = (self.yAtMove - self.yAtPress)*1
                self.moleculeRotate(0.3,yrotate,xrotate,0.0)
        elif buttons & Qt.MidButton:
            debugger.print('Mouse event - mid button')
            xshift =  0.02 * (self.xAtMove - self.xAtPress)
            yshift = -0.02 * (self.yAtMove - self.yAtPress)
            self.translate(xshift, yshift)
            self.update()
        debugger.print('Mouse event - xy', self.xAtPress,self.yAtPress)
        self.xAtPress = self.xAtMove
        self.yAtPress = self.yAtMove

    def stopAnimation(self):
        debugger.print('stopAnimation')
        if self.timer is not None:
            self.timer.stop()

    def startAnimation(self):
        debugger.print('startAnimation')
        if self.timer is not None:
            self.timer.stop()
        else:
            self.timer = QTimer(self)
            self.timer.timeout.connect(self.timeoutHandler)
            self.timer.setInterval(self.timer_interval)
        self.timer.start()

    def paintGL(self):
        debugger.print('paintGL')
        glMatrixMode(GL_MODELVIEW)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.background_colour = np.array(self.viewerTab.settings['Background colour'])/255.0
        glClearColor(*self.background_colour)
        glPushMatrix()
        glTranslatef(-self.rotation_centre[0],-self.rotation_centre[1],-self.rotation_centre[2] )
        self.drawSpheres()
        self.drawCylinders()
        if self.show_arrows:
            self.drawArrows()
        glPopMatrix()

    def drawSpheres(self):
        debugger.print('drawSpheres')
        if len(self.spheres) == 0:
            return
        for sphere in self.spheres[self.current_phase]:
            col       = sphere.colour
            rad       = sphere.radius
            (x, y, z) = sphere.position
            diffMatColour =  col * self.diffuseMaterialFactor
            ambMatColour =   col * self.ambientMaterialFactor
            specMatColour =  self.white * self.specularLightFactor
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffMatColour)
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambMatColour)
            glMaterialfv(GL_FRONT, GL_SPECULAR, specMatColour)
            glMaterialf( GL_FRONT, GL_SHININESS, self.glintMaterialFactor)
            glPushMatrix()
            glTranslatef( x, y, z )
            gluSphere(self.quadric, rad, self.sphere_slices, self.sphere_stacks)
            glPopMatrix()

    def drawCylinders(self):
        debugger.print('drawCylinders')
        if len(self.cylinders) == 0:
            return
        for cylinder in self.cylinders[self.current_phase]:
            col    = cylinder.colour
            rad    = cylinder.radius
            pos1   = cylinder.position1
            pos2   = cylinder.position2
            length = cylinder.height
            angle  = cylinder.angle
            rot    = cylinder.rotation
            diffMatColour =  col * self.diffuseMaterialFactor
            ambMatColour =   col * self.ambientMaterialFactor
            specMatColour =  self.white * self.specularLightFactor
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffMatColour)
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambMatColour)
            glMaterialfv(GL_FRONT, GL_SPECULAR, specMatColour)
            glMaterialf( GL_FRONT, GL_SHININESS, self.glintMaterialFactor)
            glPushMatrix()
            glTranslated( pos2[0], pos2[1], pos2[2] )
            glRotated(angle, rot[0], rot[1], rot[2])
            gluCylinder(self.quadric, rad, rad, length, self.cylinder_slices, self.cylinder_stacks)
            glPopMatrix()

    def drawArrows(self):
        debugger.print('drawArrows')
        if len(self.arrows) == 0:
            return
        self.arrow_colour = np.array(self.viewerTab.settings['Arrow colour'])/255.0
        for arrow,sphere in zip(self.arrows,self.spheres[self.current_phase]):
            pos       = sphere.position
            direction = arrow.direction
            length    = arrow.height
            angle     = arrow.angle
            rot       = arrow.rotation
            colour    = arrow.colour
            radius    = arrow.radius
            diffMatColour =  colour * self.diffuseMaterialFactor
            ambMatColour =   colour * self.ambientMaterialFactor
            specMatColour =  self.white * self.specularLightFactor
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffMatColour)
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambMatColour)
            glMaterialfv(GL_FRONT, GL_SPECULAR, specMatColour)
            glMaterialf( GL_FRONT, GL_SHININESS, self.glintMaterialFactor)
            glPushMatrix()
            glTranslated( pos[0], pos[1], pos[2] )
            glRotated(angle, rot[0], rot[1], rot[2])
            gluCylinder(self.quadric, radius, radius, length, self.cylinder_slices, self.cylinder_stacks)
            glTranslated( 0.0, 0.0, length )
            length = 0.1
            gluCylinder(self.quadric, 2.0*radius, 0.1*radius, length, self.cylinder_slices, self.cylinder_stacks)
            glPopMatrix()


    def resizeGL(self,w,h):
        debugger.print('resizeGL',w,h)
        self.my_width = w
        self.my_height =h
        # set projection matrix
        self.setProjectionMatrix()

    def initializeGL(self):
        debugger.print('initializeGL')
        self.quadric  = gluNewQuadric()
        gluQuadricDrawStyle(self.quadric, GLU_FILL)
        gluQuadricNormals(self.quadric, GLU_SMOOTH)
        glClearDepth(1.0)
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_NORMALIZE)
        #glDisable(GL_DITHER)
        glShadeModel(GL_SMOOTH)
        glEnable(GL_POLYGON_SMOOTH)
        glEnable(GL_POINT_SMOOTH)
        glEnable(GL_LINE_SMOOTH)
        glEnable(GL_MULTISAMPLE)
        glCullFace(GL_BACK)
        glEnable(GL_CULL_FACE)
        #glBlendFunc(GL_SRC_ALPHA_SATURATE, GL_ONE)
        glEnable(GL_BLEND)
        self.defineLights()
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.background_colour = np.array(self.viewerTab.settings['Background colour'])/255.0
        glClearColor(*self.background_colour)

    def setImageSize(self):
        maxsize = 0.0
        for sphere in self.spheres[self.current_phase]:
            pos = sphere.position
            vec = pos - self.rotation_centre
            dist = math.sqrt(np.dot(vec,vec))
            if dist > maxsize:
                maxsize = dist
        for cylinder in self.cylinders[self.current_phase]:
            pos1 = cylinder.position1
            pos2 = cylinder.position2
            vec1 = pos1 - self.rotation_centre
            dist1 = math.sqrt(np.dot(vec,vec))
            vec2 = pos2 - self.rotation_centre
            dist2 = math.sqrt(np.dot(vec,vec))
            dist = max(dist1,dist2)
            if dist > maxsize:
                maxsize = dist
        self.image_size = maxsize
        self.setProjectionMatrix()
        debugger.print('setImageSize',self.image_size)

    def setProjectionMatrix(self):
        if self.image_size is None or self.my_width is None or self.my_height is None:
            return
        self.myMakeCurrent()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        orthox = 1.1 * self.image_size * self.my_width  / min(self.my_width,self.my_height)
        orthoy = 1.1 * self.image_size * self.my_height / min(self.my_width,self.my_height)
        orthoz = 1.1 * self.image_size * max(self.my_width,self.my_height)
        debugger.print('projection', orthox, orthoy, orthoz)
        glOrtho(-orthox, orthox, -orthoy, orthoy, -orthoz, orthoz)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.matrix =  np.eye( 4, dtype=np.float32)
        # reset the current phase to the centre of the phases
        self.current_phase = int(self.number_of_phases / 2)
        debugger.print('set projection matrix ortho', orthox, orthoy, orthoz)
        debugger.print('set projection matrix image_size', self.image_size)
        debugger.print('set projection matrix current_phase', self.current_phase)

    def defineLights(self):
        debugger.print('Define Lights')
        self.myMakeCurrent()
        if self.light_switches is None:
            self.light_switches = self.viewerTab.light_switches
            glEnable(GL_LIGHTING)
            for position,lightOn,light in zip(self.light_positions,self.light_switches,self.lights):
                glLightfv(light, GL_AMBIENT,  [ self.ambientLightFactor,  self.ambientLightFactor,  self.ambientLightFactor,  1.0 ] )
                glLightfv(light, GL_DIFFUSE,  [ self.diffuseLightFactor,  self.diffuseLightFactor,  self.diffuseLightFactor,  1.0 ] )
                glLightfv(light, GL_SPECULAR, [ self.specularLightFactor, self.specularLightFactor, self.specularLightFactor, 1.0 ] )
                glLightfv(light, GL_POSITION, position )
                if self.linearAttenuation:
                    glLight(light, GL_LINEAR_ATTENUATION, 1.0)
        for lightOn, light in zip(self.light_switches,self.lights):
            if lightOn:
                glEnable(light)
            else:
                glDisable(light)

    def setRotationCentre(self, pos):
        debugger.print('set rotation centre',pos)
        self.rotation_centre = np.array(pos)

    def createArrays(self, nphases):
        debugger.print('createArrays')
        #  Create an empty list for each phase
        # self.spheres    = [ [] for i in range(nphases) ]
        # self.cylinders  = [ [] for i in range(nphases) ]
        #self.spheres    = deque( deque() for i in range(nphases) )
        #self.cylinders  = deque( deque() for i in range(nphases) )
        self.spheres.clear()
        self.cylinders.clear()
        for i in range(nphases):
            self.spheres.append(deque())
            self.cylinders.append(deque())
        self.number_of_phases = nphases
        self.current_phase = int(self.number_of_phases / 2)

    def deleteSpheres(self):
        debugger.print('deleteSpheres')
        self.spheres.clear()

    def deleteCylinders(self):
        debugger.print('deleteCylinders')
        self.cylinders.clear()

    def deleteArrows(self):
        debugger.print('deleteArrows')
        #self.arrows = []
        self.arrows.clear()

    def addArrows(self, colour, radius, direction, length, phase=0):
        # There is no phase requirement for arrows - they are just displacements
        self.arrows.append( Arrow(colour, radius, direction, length) )

    def addCylinder(self, colour, radius, pos1, pos2, phase=0):
        self.cylinders[phase].append( Cylinder(colour,radius,pos1,pos2) )

    def addSphere(self, colour, radius, pos, phase=0):
        self.spheres[phase].append( Sphere(colour,radius,pos) )

class Sphere():
    def __init__(self, colour, radius, position):
        self.colour= np.array(colour)/255.0
        self.radius = radius
        self.position = np.array(position)

class Cylinder():
    def __init__(self, colour, radius, pos1, pos2):
        a = np.array(pos1)
        b = np.array(pos2)
        self.colour= np.array(colour)/255.0
        self.radius = radius
        self.position1 = a
        self.position2 = b
        z = np.array( [0.0, 0.0, 1.0] )
        p = a - b
        t = np.cross(z, p)
        height = math.sqrt(np.dot(p,p))
        angle = 180.0 / PI * math.acos( np.dot(z,p) / height )
        self.angle = angle
        self.height = height
        self.rotation = t

class Arrow():
    def __init__(self, colour, radius, direction, scale):
        # Adds a displacement to each atom
        self.colour= np.array(colour)/255.0
        p = np.array(direction)
        z = np.array( [0.0, 0.0, 1.0] )
        t = np.cross(z, p)
        height = math.sqrt(np.dot(p,p)) + 1.0E-8
        angle = 180.0 / PI * math.acos( np.dot(z,p) / height )
        p = p / height
        self.angle = angle
        self.height = height*scale
        self.rotation = t
        self.direction = p
        self.radius = radius

