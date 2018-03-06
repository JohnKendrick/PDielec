import numpy as np
import math
from OpenGL.GL       import *
from OpenGL.GLU      import *
from OpenGL.GLUT     import *
from PyQt5.QtWidgets import QOpenGLWidget
from PyQt5.QtCore    import Qt
from Python.Constants import PI

class OpenGLWidget(QOpenGLWidget):

    def __init__(self, parent):
        QOpenGLWidget.__init__(self, parent)
        self.setMinimumSize(640, 680)
        self.lightingOn = True
        self.diffuseMaterialFactor  = 0.9
        self.ambientMaterialFactor  = 0.4
        self.specularMaterialFactor = 0.7
        self.glintMaterialFactor    = 100.0
        self.diffuseLightFactor     = 0.8
        self.ambientLightFactor     = 0.6
        self.specularLightFactor    = 1.0
        self.linearAttenuation      = True
        self.lightSwitches          = [False]*8
        self.lightSwitches[0]       = True
        self.lightSwitches[1]       = True
        self.image_size             = 15.0
        self.white                  = np.array( [1.0, 1.0, 1.0, 1.0] )
        self.sphere_colours         = []
        self.sphere_radii           = []
        self.sphere_positions       = []
        self.sphere_slices          = 15
        self.sphere_stacks          = 15
        self.cylinder_colours       = []
        self.cylinder_radii         = []
        self.cylinder_pos1          = []
        self.cylinder_pos2          = []
        self.cylinder_slices        = 10
        self.cylinder_stacks        =  2
        self.cylinder_angles        = []
        self.cylinder_heights       = []
        self.cylinder_rotations     = []
        self.arrow_directions       = []
        self.arrow_angles           = []
        self.arrow_heights          = []
        self.arrow_rotations        = []
        self.arrow_scaling          = 10
        self.arrow_colour           = np.array( [ 0.0, 1.0, 0.0 ] )
        self.arrow_radius           = 0.05
        
        self.show_full_screen       = False
        self.zoom_factor            = 1.0
        self.rotation_centre         = np.array( [0.0, 0.0, 0.0] )
        self.quadric                = gluNewQuadric()
        gluQuadricDrawStyle(self.quadric, GLU_FILL)
        gluQuadricNormals(self.quadric, GLU_SMOOTH)
        self.setFocusPolicy(Qt.StrongFocus)
        self.matrix =  np.eye( 4, dtype=np.float32)
#        format = QSurfaceFormat()
#        format.setVersion(2,0)
#        format.setBufferSize(24)
#        format.setStencilBufferSize(8)
#        QSurfaceFormat.sefDefaultFormat(format)

    def moleculeRotate(self,scale,x,y,z):
        # Rotate molecular frame using glortho
        # Create a 4x4 matrix for the model view
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glRotatef(scale*x, 1.0, 0.0, 0.0)
        glRotatef(scale*y, 0.0, 1.0, 0.0)
        glMultMatrixf(self.matrix)
        glGetFloatv(GL_MODELVIEW_MATRIX, self.matrix)
        self.update()

    def debug_trace(self):
        from PyQt5.QtCore import pyqtRemoveInputHook
        import pdb
        pyqtRemoveInputHook()
        pdb.set_trace()

    def keyPressEvent(self, event):
        key = event.key()
        modifiers = event.modifiers()
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
        elif key == Qt.Key_F:
            self.show_full_screen = not self.show_full_screen
            print('key is f', self.show_full_screen)
            if self.show_full_screen:
                self.showFullScreen()
            else:
                self.showNormal()
        elif key == Qt.Key_Home:
            print('home key pressed')
            glMatrixMode(GL_MODELVIEW)
            glLoadIdentity()
            self.matrix =  np.eye( 4, dtype=np.float32)
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
        if zoom > 0:
            self.zoom_factor *= 1.03
        else:
            self.zoom_factor *= 0.97
        self.update()

    def mouseMoveEvent(self, event):
        self.xAtMove = event.x()
        self.yAtMove = event.y()
        modifiers = event.modifiers()
        if modifiers & Qt.ShiftModifier or modifiers & Qt.ControlModifier:
            # handle zoom
            xzoom = +1.0*(self.xAtMove - self.xAtPress)
            yzoom = -1.0*(self.yAtMove - self.yAtPress)
            self.zoom(xzoom+yzoom)
        else:
            # handle ordinary rotation
            xrotate = (self.xAtMove - self.xAtPress)*1
            yrotate = (self.yAtMove - self.yAtPress)*1
            self.moleculeRotate(1.0,yrotate,xrotate,0.0)
        self.xAtPress = self.xAtMove
        self.yAtPress = self.yAtMove

    def paintGL(self):
        glMatrixMode(GL_MODELVIEW)
        # For some reason the paintGL function seems to reset the modelview matrix
        glLoadMatrixf(self.matrix)
        # For the same reason we need to apply the zoom factor here
        glScalef(self.zoom_factor, self.zoom_factor, self.zoom_factor)
        # Clear the screen
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glPushMatrix()
        glTranslatef(-self.rotation_centre[0],-self.rotation_centre[1],-self.rotation_centre[2] )
        self.drawSpheres()
        self.drawCylinders()
        self.drawArrows()
        glPopMatrix()

    def drawSpheres(self):
        for (col,rad, (x, y, z)) in zip(self.sphere_colours, self.sphere_radii, self.sphere_positions ):
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
        for (col,rad, pos1, pos2, length, angle, rot) in zip(self.cylinder_colours, self.cylinder_radii, self.cylinder_pos1, self.cylinder_pos2, self.cylinder_heights, self.cylinder_angles, self.cylinder_rotations):
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
        for ( pos, direction, length, angle, rot) in zip( self.sphere_positions, self.arrow_directions, self.arrow_heights, self.arrow_angles, self.arrow_rotations):
            diffMatColour =  self.arrow_colour * self.diffuseMaterialFactor
            ambMatColour =   self.arrow_colour * self.ambientMaterialFactor
            specMatColour =  self.white * self.specularLightFactor
            glMaterialfv(GL_FRONT, GL_DIFFUSE, diffMatColour)
            glMaterialfv(GL_FRONT, GL_AMBIENT, ambMatColour)
            glMaterialfv(GL_FRONT, GL_SPECULAR, specMatColour)
            glMaterialf( GL_FRONT, GL_SHININESS, self.glintMaterialFactor)
            length *= self.arrow_scaling
            glPushMatrix()
            glTranslated( pos[0], pos[1], pos[2] )
            glRotated(angle, rot[0], rot[1], rot[2])
            gluCylinder(self.quadric, self.arrow_radius, self.arrow_radius, length, self.cylinder_slices, self.cylinder_stacks)
            glTranslated( 0.0, 0.0, length )
            length = 0.1
            gluCylinder(self.quadric, 2.0*self.arrow_radius, 0.1*self.arrow_radius, length, self.cylinder_slices, self.cylinder_stacks)
            glPopMatrix()


    def resizeGL(self,w,h):
        print('resizeGL',w,h)
        self.width = w
        self.height =h
        # set projection matrix
        self.setProjectionMatrix()

    def initializeGL(self):
        print('Initialise')
        self.makeCurrent()
        glutInit('3D Visualiser')
        glClearDepth(1.0)              
        glDepthFunc(GL_LESS)
        glEnable(GL_DEPTH_TEST)
        glEnable(GL_NORMALIZE)
        glDisable(GL_DITHER)
        glShadeModel(GL_SMOOTH)
        self.defineLights()

    def setImageSize(self):
        maxsize = 0.0
        for sphere in self.sphere_positions:
            vec = sphere - self.rotation_centre
            dist = math.sqrt(np.dot(vec,vec))
            if dist > maxsize:
                maxsize = dist
        for pos in self.cylinder_pos1:
            vec = pos - self.rotation_centre
            dist = math.sqrt(np.dot(vec,vec))
            if dist > maxsize:
                maxsize = dist
        for pos in self.cylinder_pos2:
            vec = pos - self.rotation_centre
            dist = math.sqrt(np.dot(vec,vec))
            if dist > maxsize:
                maxsize = dist
        print('maximum  size', maxsize)
        self.image_size = maxsize
        self.setProjectionMatrix()
        

    def setProjectionMatrix(self):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()                    
        orthox = 1.1 * self.image_size * self.width  / max(self.width,self.imageSize)
        orthoy = 1.1 * self.imageS_size * self.height / max(self.width,self.imageSize)
        orthoz = 1.1 * self.image_size
        glOrtho(-orthox, orthox, -orthoy, orthoy, -orthoz, orthoz)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        self.matrix =  np.eye( 4, dtype=np.float32)
        self.zoom_factor = 1.0

    def defineLights(self):
        glEnable(GL_LIGHTING)
        d = 200.0
        positions = [ [-d, d, d], [ d,-d, d], [-d,-d, d], [ d, d, d], [-d, d,-d], [ d,-d,-d], [-d,-d,-d], [ d, d,-d] ]
        lights = [ GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3, GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7 ]
        for position,lightOn,light in zip(positions,self.lightSwitches,lights):
            glLightfv(light, GL_AMBIENT,  [ self.ambientLightFactor,  self.ambientLightFactor,  self.ambientLightFactor,  1.0 ] )
            glLightfv(light, GL_DIFFUSE,  [ self.diffuseLightFactor,  self.diffuseLightFactor,  self.diffuseLightFactor,  1.0 ] )
            glLightfv(light, GL_SPECULAR, [ self.specularLightFactor, self.specularLightFactor, self.specularLightFactor, 1.0 ] )
            glLightfv(light, GL_POSITION, position )
            if self.linearAttenuation:
                glLight(light, GL_LINEAR_ATTENUATION, 1.0)
            if lightOn:
                glEnable(light)
            else:
                glDisable(light)

    def setRotationCentre(self, pos):
        self.rotation_centre = np.array(pos)

    def deleteSpheres(self):
        self.sphere_colours = []
        self.sphere_radii = []
        self.sphere_positions = []

    def deleteCylinders(self):
        self.cylinder_colours = []
        self.cylinder_radii = []
        self.cylinder_pos1 = []
        self.cylinder_pos2 = []
        self.cylinder_angles = []
        self.cylinder_heights = []
        self.cylinder_rotations = []

    def deleteDisplacements(self):
        # This is called subsequently to all arrows having been added, it sets all the arrow directions for every atom
        self.arrow_directions = []
        self.arrow_angles = []
        self.arrow_heights = []
        self.arrow_rotations = []

    def addDisplacement(self, direction):
        # Adds a displacement to each atom
        p = np.array(direction)
        z = np.array( [0.0, 0.0, 1.0] )
        t = np.cross(z, p)
        height = math.sqrt(np.dot(p,p)) + 1.0E-8
        angle = 180.0 / PI * math.acos( np.dot(z,p) / height )
        p = p / height
        self.arrow_angles.append(angle)
        self.arrow_heights.append(height)
        self.arrow_rotations.append(t)
        self.arrow_directions.append(p)

    def addCylinder(self, colour, radius, pos1, pos2):
        a = np.array(pos1)
        b = np.array(pos2)
        self.cylinder_colours.append(np.array(colour))
        self.cylinder_radii.append(radius)
        self.cylinder_pos1.append(a)
        self.cylinder_pos2.append(b)
        z = np.array( [0.0, 0.0, 1.0] )
        p = a - b
        t = np.cross(z, p)
        height = math.sqrt(np.dot(p,p))
        angle = 180.0 / PI * math.acos( np.dot(z,p) / height )
        self.cylinder_angles.append(angle)
        self.cylinder_heights.append(height)
        self.cylinder_rotations.append(t)

    def addSphere(self, colour, radius, pos):
        self.sphere_colours.append(np.array(colour))
        self.sphere_radii.append(radius)
        self.sphere_positions.append(np.array(pos))
