#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""OpenGLWidget module."""

import math
import os
from collections import deque

import numpy as np
import OpenGL
from OpenGL.GL import (
    GL_ALL_ATTRIB_BITS,
    GL_AMBIENT,
    GL_BACK,
    GL_BLEND,
    GL_COLOR_BUFFER_BIT,
    GL_CULL_FACE,
    GL_DEPTH_BUFFER_BIT,
    GL_DEPTH_TEST,
    GL_DIFFUSE,
    GL_FRONT,
    GL_LESS,
    GL_LIGHT0,
    GL_LIGHT1,
    GL_LIGHT2,
    GL_LIGHT3,
    GL_LIGHT4,
    GL_LIGHT5,
    GL_LIGHT6,
    GL_LIGHT7,
    GL_LIGHTING,
    GL_LINE_SMOOTH,
    GL_LINEAR_ATTENUATION,
    GL_MODELVIEW,
    GL_MODELVIEW_MATRIX,
    GL_MULTISAMPLE,
    GL_NO_ERROR,
    GL_NORMALIZE,
    GL_POINT_SMOOTH,
    GL_POLYGON_SMOOTH,
    GL_POSITION,
    GL_PROJECTION,
    GL_PROJECTION_MATRIX,
    GL_SHININESS,
    GL_SMOOTH,
    GL_SPECULAR,
    GL_VIEWPORT,
    glClear,
    glClearColor,
    glClearDepth,
    glCullFace,
    glDepthFunc,
    glDisable,
    glEnable,
    glGetDoublev,
    glGetError,
    glGetFloatv,
    glGetIntegerv,
    glLight,
    glLightfv,
    glLoadIdentity,
    glMaterialf,
    glMaterialfv,
    glMatrixMode,
    glMultMatrixf,
    glOrtho,
    glPopAttrib,
    glPopMatrix,
    glPushAttrib,
    glPushMatrix,
    glRotated,
    glRotatef,
    glScalef,
    glShadeModel,
    glTranslated,
    glTranslatef,
)
from OpenGL.GLU import (
    GLU_FILL,
    GLU_SMOOTH,
    gluCylinder,
    gluLookAt,
    gluNewQuadric,
    gluProject,
    gluQuadricDrawStyle,
    gluQuadricNormals,
    gluSphere,
)
from qtpy.QtCore import Qt, QTimer
from qtpy.QtGui import QFont, QPainter, QSurfaceFormat
from qtpy.QtWidgets import QMessageBox, QOpenGLWidget

from PDielec.Utilities import Debug

# The following lines seem to fix a problem when running on low end machines
OpenGL.USE_ACCELERATE = False
OpenGL.ERROR_CHECKING = False
# end of low-end machine fix


class OpenGLWidget(QOpenGLWidget):
    """A widget for the Qt framework specialized in OpenGL graphics rendering.

    This widget is designed to act as a capable 3D rendering environment
    for OpenGL within a Qt application. It supports a wide range of 3D
    visualization tasks, including rendering spheres, cylinders, and arrows
    to represent molecular structures or other 3D objects.

    Parameters
    ----------
    parent : QWidget
        The parent widget.
    debug : bool, optional
        Enables debug output if True.

    Attributes
    ----------
    format : QSurfaceFormat
        Specifies the format options for the rendering context. Common settings include
        samples for multisampling and the swap behavior for double buffering.
    spheres : deque
        A queue that holds sphere objects to be rendered. These can represent atoms or other spherical entities.
    cylinders : deque
        A queue that holds cylinder objects representing bonds or connections between spheres.
    texts : deque
        A queue that holds text objects.
    arrows : deque
        A queue that holds arrow objects, useful for indicating directions or forces.
    light_positions : list
        Defines positions for up to eight lights in the 3D space.
    lights : list
        References to OpenGL light objects, typically provided as constants like GL_LIGHT0, GL_LIGHT1, etc.

    Methods
    -------
    enterEvent(event)
        Handles the mouse entering the widget's area.
    showArrows(show)
        Toggles the visibility of arrows.
    timeoutHandler()
        Advances the current animation phase or shifts direction at bounds.
    myMakeCurrent()
        Ensures the OpenGL context is current and checks for errors.
    moleculeRotate(scale, x, y, z)
        Applies a rotation to the molecule model.
    keyPressEvent(event)
        Custom key press event handler for zooming and rotating the model.
    save_movie(filename)
        Saves a rendered animation to a movie file.
    snapshot(filename)
        Takes a snapshot of the current OpenGL frame and saves it as an image.
    translate(x, y)
        Translates the scene.
    wheelEvent(event)
        Handles mouse wheel events for zooming.
    mousePressEvent(event), mouseReleaseEvent(event), mouseMoveEvent(event)
        Handles mouse interactions for rotating and translating the scene.
    zoom(zoom)
        Applies zoom to the scene based on input.
    stopAnimation()
        Stops the animation timer.
    startAnimation()
        Starts or restarts the animation timer.
    paintGL()
        The main rendering function called by the Qt framework to draw the OpenGL scene.
    drawSpheres(), drawCylinders(), drawArrows(), drawTexts()
        Functions responsible for drawing the respective geometric objects.
    resizeGL(w, h)
        Adjusts the viewport and projection matrix on widget resize.
    initializeGL()
        Initializes OpenGL states and geometry.
    setImageSize()
        Calculates the ideal image size for the viewport.
    setProjectionMatrix()
        Sets up the projection matrix based on the current image size and viewport dimensions.
    defineLights()
        Configures the light sources.
    setRotationCentre(pos)
        Sets the center of rotation for the scene.
    createArrays(nphases)
        Initializes storage for drawable objects across animation phases.
    deleteSpheres(), deleteCylinders(), deleteArrows(), deleteTexts()
        Convenience methods for clearing the drawable objects.
    addArrows(colour, radius, direction, length, phase=0)
        Adds an arrow to the scene.
    addCylinder(colour, radius, pos1, pos2, phase=0)
        Adds a cylinder to the scene.
    addSphere(colour, radius, pos, phase=0)
        Adds a sphere to the scene.
    addText(colour, size, pos, phase=0)
        Adds a text to the scene.

    """

    format = QSurfaceFormat()
    #format.setDepthBufferSize(24)
    #format.setStencilBufferSize(8)
    format.setSamples(6)
    format.setSwapBehavior(QSurfaceFormat.DoubleBuffer)
    QSurfaceFormat.setDefaultFormat(format)

    def __init__(self, parent, debug=False):
        """Initialise the OpenGLWidget object.

        This method initializes an OpenGL Widget with several default properties such as lighting, material factors, and visualization elements like spheres, cylinders, and arrows. It also sets up the widget's size and various attributes related to phase animation and visualization settings.

        Parameters
        ----------
        parent : QWidget
            The parent QWidget for this OpenGLWidget.
        debug : bool, optional
            A boolean flag to set debugging mode. The default is False.

        Notes
        -----
        The constructor initializes various lists and attributes to default values. It also inherits from QOpenGLWidget and thus requires calling the superclass constructor.
        Attributes like lighting configurations, material factors (diffuse, ambient, specular), phase control attributes (current phase, phase direction, number of phases), and geometric object counts and properties (spheres, cylinders, arrows) are initialized.
        Moreover, the widget's minimum size is set, and default values for lighting intensity factors, color, phase animation, geometrical properties, light positions, and switches are provided.
        Attributes initialized include lighting and material factors, object lists for spheres, cylinders, and arrows to visualize, phase control, lighting positions, canvas size, and orientation, among others.

        Example
        -------
        To create an instance of the OpenGLWidget with debugging enabled, you could use:

        >>> parent_widget = SomeQWidget()
        >>> opengl_widget = OpenGLWidget(parent=parent_widget, debug=True)

        """        
        QOpenGLWidget.__init__(self, parent)
        global debugger
        debugger = Debug(debug,"OpenGLWidget")
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
        self.texts                  = deque()
        self.current_phase          = 0
        self.phase_direction        = 1
        self.number_of_phases       = 1
        self.sphere_slices          = 20
        self.sphere_stacks          = 20
        self.cylinder_slices        = 8
        self.cylinder_stacks        = 2
        self.timer                  = None
        self.show_arrows            = True
        self.show_orientation       = False
        self.timer_interval         = 60
        self.my_width               = None
        self.my_height              = None
        self.background_colour      = None
        self.show_full_screen       = False
        self.writer                 = None
        self.rotation_centre        = np.array( [0.0, 0.0, 0.0] )
        self.matrix                 =  np.eye( 4, dtype=np.float32)
        # The orientation definitions are up, across, out vectors respectively
        self.orientation_definitions = { "x": ( [0,1,0], [0,0,1], [1,0,0] ),
                                         "y": ( [0,0,1], [1,0,0], [0,1,0] ),
                                         "z": ( [1,0,0], [0,1,0], [0,0,1] ),
                                         "X": ( [1,0,0], [0,1,0], [0,0,1] ),
                                         "Y": ( [0,1,0], [0,0,1], [1,0,0] ),
                                         "Z": ( [0,0,1], [1,0,0], [0,1,0] ),
        }
        self.orientation             = "z"
        self.current_orientation     = self.orientation_definitions[self.orientation]
        self.hkl                     = (1,0,0)
        self.cell                    = None
        self.light_switches          = None
        d = 200.0
        self.light_positions = [ [-d, d, d], [ d,-d, d], [-d,-d, d], [ d, d, d], [-d, d,-d], [ d,-d,-d], [-d,-d,-d], [ d, d,-d] ]
        self.lights = [ GL_LIGHT0, GL_LIGHT1, GL_LIGHT2, GL_LIGHT3, GL_LIGHT4, GL_LIGHT5, GL_LIGHT6, GL_LIGHT7 ]

    def enterEvent(self, event):
        """Handle the entering of an event.

        Parameters
        ----------
        event : QEvent
            The event that triggers the enterEvent.

        Returns
        -------
        None

        Notes
        -----
        This function primarily ensures that the focus is set accordingly when an event occurs, and then it passes the event on to the viewerTab's own enterEvent handler for further processing.

        """        
        debugger.print("enter event")
        self.setFocus()
        self.viewerTab.enterEvent(event)

    def showArrows(self, show):
        """Set the visibility of arrows.

        Parameters
        ----------
        show : bool
            A boolean indicating whether to show arrows (True) or not (False).

        Returns
        -------
        None

        Notes
        -----
        This function updates the `show_arrows` attribute based on the "show" parameter and prints the current status to the debugger.

        """        
        debugger.print("show arrows",show)
        self.show_arrows = show

    def timeoutHandler(self):
        """Handle the timeout event in a cyclic phase system by updating the current phase and direction.

        This function updates the `current_phase` by advancing or reversing it based on the `phase_direction`.
        When the end or the start of the phases is reached, the direction is reversed, ensuring cyclic behavior
        through the phases. After updating the phase and direction, the function logs the current phase and
        updates the system state.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - `self.current_phase` is the index of the current phase in the cycle.
        - `self.number_of_phases` indicates the total number of available phases.
        - `self.phase_direction` is either +1 or -1, indicating the direction of phase progression.
        - The `update` method is called at the end to apply any changes made to the state as a result of this timeout event.

        """        
        self.current_phase += self.phase_direction
        if self.current_phase >= self.number_of_phases:
            self.phase_direction = -1
            self.current_phase = self.number_of_phases -2
        elif self.current_phase < 0:
            self.current_phase = 1
            self.phase_direction = +1
        debugger.print("Timeout - phase", self.current_phase)
        self.update()

    def myMakeCurrent(self):
        """Make the current OpenGL context active and clear any errors.

        This method activates the current OpenGL context by calling `makeCurrent`,
        then repeatedly calls `glGetError` to clear any accumulated OpenGL errors.

        Notes
        -----
        This function requires that an OpenGL context is available and that `makeCurrent`
        and `glGetError` are correctly implemented in the class this method belongs to.
        `makeCurrent` should make the OpenGL context current, and `glGetError` should
        return the current OpenGL error code, with `GL_NO_ERROR` indicating no error.

        This method does not return any value or raise exceptions. It is designed to
        ensure that the OpenGL context is error-free after making it current.

        """        
        self.makeCurrent()
        err = glGetError()
        while err != GL_NO_ERROR:
            err = glGetError()

    def moleculeRotate(self,scale,x,y,z):
        """Rotate a molecule object around the specified axes.

        This method applies a rotation transformation to a molecule object
        around the x, y, and z axes, scaled by a given factor.

        Parameters
        ----------
        scale : float
            The scaling factor for the rotation. Represents how much the rotation is amplified.
        x : float
            The degree of rotation around the x-axis (degrees).
        y : float
            The degree of rotation around the y-axis (degrees).
        z : float
            Degree of rotation around the z-axis (degrees).

        Returns
        -------
        None

        Notes
        -----
        - This function assumes that the GL_MODELVIEW matrix mode is active and manipulates the current matrix stack accordingly.
        - It first saves the current model view matrix, then resets it to the identity.
        - Rotations for the x and y axes are applied in sequence, scaled by the provided `scale` factor.
        - Finally, the original matrix is restored and multiplied by the new rotation transformations.
        - The function uses OpenGL calls like `glRotatef` and `glMultMatrixf`, making sure the graphical representation of the molecule is updated to reflect the rotation.

        Examples
        --------
        >>> mol_object = Molecule()
        >>> mol_object.moleculeRotate(45, 1, 0, 0)
        This would rotate the molecule 45 degrees around the x-axis.

        See Also
        --------
        glRotatef : Rotates the current matrix a specified number of degrees about a specified axis.
        glMultMatrixf : Multiplies the current matrix by the specified matrix.
        glLoadIdentity : Replaces the current matrix with the identity matrix.

        Raises
        ------
        This function does not explicitly raise exceptions but relies on the OpenGL context and the `update` method's implementation.

        """        
        # Rotate molecular frame using glortho
        # Create a 4x4 matrix for the model view
        debugger.print("molecule rotation",x,y,z)
        self.myMakeCurrent()
        glMatrixMode(GL_MODELVIEW)
        self.matrix=glGetFloatv(GL_MODELVIEW_MATRIX)
        old_matrix = self.matrix[:3,:3]
        up,across,out = self.current_orientation
        debugger.print("up",up)
        debugger.print("across",across)
        debugger.print("out",out)
        glRotatef(scale*x, across[0], across[1], across[2])
        glRotatef(scale*y,     up[0],     up[1],     up[2])
        glRotatef(scale*z,    out[0],    out[1],    out[2])
        self.matrix=glGetFloatv(GL_MODELVIEW_MATRIX)
        new_matrix = self.matrix[:3,:3]
        new_up     = np.dot(new_matrix,np.dot(old_matrix.T,up))
        new_across = np.dot(new_matrix,np.dot(old_matrix.T,across))
        new_out    = np.dot(new_matrix,np.dot(old_matrix.T,out))
        new_up     = new_up / np.linalg.norm(new_up)
        new_across = new_across / np.linalg.norm(new_across)
        new_out    = new_out / np.linalg.norm(new_out)
        self.current_orientation = (new_up, new_across, new_out)
        debugger.print("new_up",new_up)
        debugger.print("new_across",new_across)
        debugger.print("new_out",new_out)
        uvw_up     = self.cell.convert_xyz_to_integer_abc(new_up)
        uvw_across = self.cell.convert_xyz_to_integer_abc(new_across)
        uvw_out    = self.cell.convert_xyz_to_integer_abc(new_out)
        self.current_uvw_orientation = (uvw_up, uvw_across, uvw_out)
        debugger.print("uvw_up",uvw_up)
        debugger.print("uvw_across",uvw_across)
        debugger.print("uvw_out",uvw_out)
        self.update()

    def show_help_dialog(self):
        """Show the help dialogue message.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        help_string = """
  Key      	Description                                         
  ------	--------------------------------------------------  
  ⇒      	: Rotate around screen vertical axis 5 degrees       
  Shift⇒ 	: Rotate around screen vertical axis  45 degrees     
  ⇐      	: Rotate around screen vertical axis -5 degrees      
  Shift⇐ 	: Rotate around screen vertical axis  -45 degrees    
  ⇑      	: Rotate around screen horizontal axis 5 degrees     
  Shift⇑ 	: Rotate around screen horizontal axis  45 degrees   
  ⇓      	: Rotate around screen horizontal axis -5 degrees    
  Shift⇓ 	: Rotate around screen horizontal axis  -45 degrees  
  l      	: Rotate cell anticlockwise  5 degrees               
  L      	: Rotate cell anticlockwise  45 degrees              
  r      	: Rotate cell clockwise  5 degrees                   
  R      	: Rotate cell clockwise  45 degrees                 
  x      	: Show the cell with the x-axis normal to the screen   
  X      	: Show the cell with the x-axis vertical on the screen          
  y      	: Show the cell with the y-axis normal to the screen   
  Y      	: Show the cell with the y-axis vertical on the screen          
  z      	: Show the cell with the z-axis normal to the screen   
  Z      	: Show the cell with the Z-axis vertical on the screen          
  a      	: Show the cell with the a-axis normal to the screen   
  A      	: Show the cell with the a-axis vertical on the screen          
  {ctrl}a	: Show the cell with the reciprocal a-axis normal to the screen  
  {ctrl}A	: Show the cell with the reciprocal a-axis vertical on the screen
  b      	: Show the cell with the b-axis normal to the screen
  B      	: Show the cell with the b-axis vertical on the screen 
  {ctrl}b	: Show the cell with the reciprocal b-axis normal to the screen
  {ctrl}B	: Show the cell with the reciprocal b-axis vertical on the screen
  c      	: Show the cell with the c-axis normal to the screen
  C      	: Show the cell with the c-axis vertical on the screen
  {ctrl}c	: Show the cell with the reciprocal c-axis normal to the screen
  {ctrl}C	: Show the cell with the reciprocal c-axis vertical on the screen
  s      	: Show (hkl) surface normal & [uvw] vertical on the screen
  S      	: Show (hkl) surface vertical & [uvw] normal on the screen
"""
        msg = QMessageBox(self)
        msg.setWindowTitle("Keyboard Shortcuts Help")
        msg.setText(help_string)
        msg.exec()
      
   

    def keyPressEvent(self, event):
        """Handle keyboard events for specific actions.

        This method is triggered whenever a keyboard event occurs. 
        It interprets key presses to perform various actions such as rotating the molecule, zooming in or out, 
        saving the current state, toggling between full screen mode, and resetting the view. 
 
        Parameters
        ----------
        event : QKeyEvent
            The event object that contains information about the key press event.

        Returns
        -------
        None

        Notes
        -----
        This function checks for specific keys and their combinations with modifiers like Shift and Control to determine the required action. 
        For rotation and zoom actions, the amount of change is determined based on whether 
        any modifiers are pressed (larger rotation if Shift or Control is held). 
        Saving, toggling the full screen mode, and resetting the view do not require modifiers.

        The following keys perform actions as described:
        - Left/Right Arrow: Rotate the molecule along the Y-axis.
        - Up/Down Arrow: Rotate the molecule along the X-axis.
        - Plus/Minus: Zoom in or out.
        - P: Save current state/configuration.
        - F: Toggle full screen mode.
        - Home: Reset the view to the initial state.

        Modifiers:
        - Shift: Increases the rotation amount.
        - Control: Increases the rotation amount.

        Attributes modified:
        - This method modifies various attributes of the instance based on the key pressed, including rotation angles, zoom levels, screen mode, and view matrix.

        Requires
        --------
        This function requires `Qt`, `GL_MODELVIEW`, `glLoadIdentity`, `np`, and the instance methods 
        `myMakeCurrent`, `moleculeRotate`, `zoom`, `save_movie`, `showFullScreen`, and `showNormal` to be properly defined. 
        It also assumes an external debugger tool for logging.

        Raises
        ------
        This method does not explicitly raise any exceptions but might propagate exceptions from methods it calls.

        """        
        self.myMakeCurrent()
        key = event.key()
        modifiers = event.modifiers()
        control = False
        shift = False
        if modifiers & Qt.ShiftModifier:
            shift = True
        if modifiers & Qt.ControlModifier:
            control = True
        debugger.print("kepressevent",key,modifiers,control,shift)
        amount = 45.0 if modifiers & Qt.ShiftModifier or modifiers & Qt.ControlModifier else 5.0
        if key == Qt.Key_Left:
            self.moleculeRotate(+amount,0.0,1.0,0.0)
        elif key == Qt.Key_Right:
            self.moleculeRotate(-amount,0.0,1.0,0.0)
        elif key == Qt.Key_Up:
            self.moleculeRotate(-amount,1.0,0.0,0.0)
        elif key == Qt.Key_Down:
            self.moleculeRotate(+amount,1.0,0.0,0.0)
        elif key == Qt.Key_R:
            self.moleculeRotate(-amount,0.0,0.0,1.0)
        elif key == Qt.Key_L:
            self.moleculeRotate(+amount,0.0,0.0,1.0)
        elif key == Qt.Key_H:
            self.show_help_dialog()
        elif key == Qt.Key_Plus:
            self.zoom(+1.0)
        elif key == Qt.Key_Minus:
            self.zoom(-1.0)
        elif key == Qt.Key_A:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("A")
                if modifiers & Qt.ControlModifier:
                    self.set_orientation("A*")
            elif modifiers & Qt.ControlModifier:
                self.set_orientation("a*")
            else:
                self.set_orientation("a")
        elif key == Qt.Key_B:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("B")
                if modifiers & Qt.ControlModifier:
                    self.set_orientation("B*")
            elif modifiers & Qt.ControlModifier:
                self.set_orientation("b*")
            else:
                self.set_orientation("b")
        elif key == Qt.Key_C:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("C")
                if modifiers & Qt.ControlModifier:
                    self.set_orientation("C*")
            elif modifiers & Qt.ControlModifier:
                self.set_orientation("c*")
            else:
                self.set_orientation("c")
        elif key == Qt.Key_X:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("X")
            else:
                self.set_orientation("x")
        elif key == Qt.Key_Y:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("Y")
            else:
                self.set_orientation("y")
        elif key == Qt.Key_Z:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("Z")
            else:
                self.set_orientation("z")
        elif key == Qt.Key_S:
            if modifiers & Qt.ShiftModifier:
                self.set_orientation("S")
            else:
                self.set_orientation("s")
        elif key == Qt.Key_P:
            self.save_movie("movie.mp4")
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
            self.set_orientation("z")
            debugger.print("Home key", self.current_phase)
            self.update()

    def set_orientation(self,orientation):
        """Set the molecular orientation.

        The orientation can be one of 'x', 'y', 'z', 'a', 'b', 'c'
        The orientation is set using gluLookAt to define the modelview matrix
        The rotation parts of the transformation matrix are maintained
        The translation part is set to unity.

                             Points out    Points up

                   "x" means     x-axis       y-axis
                   "y" means     y-axis       z-axis
                   "z" means     z-axis       x-axis
                   "a" means     a-axis       aXb
                   "b" means     b-axis       bXc
                   "c" means     c-axis       cXa

        The use of a capital letter means the axes are reversed

        Parameters
        ----------
        orientation : str
            The required orientation

        Returns
        -------
        None

        """        
        debugger.print("set_orientation", orientation)
        up,across,out = self.orientation_definitions[orientation]
        debugger.print("out   :",out)
        debugger.print("across:",across)
        debugger.print("up    :",up)
        self.orientation = orientation
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glPushMatrix()
        gluLookAt(200*out[0], 200*out[1], 200*out[2], 0, 0, 0, up[0], up[1], up[2])
        self.matrix = glGetFloatv(GL_MODELVIEW_MATRIX)
        glPopMatrix()
        self.matrix[0,3] = 0
        self.matrix[1,3] = 0
        self.matrix[2,3] = 0
        self.matrix[3,0] = 0
        self.matrix[3,1] = 0
        self.matrix[3,2] = 0
        self.matrix[3,3] = 1
        glMultMatrixf(self.matrix)
        uvw_out    = self.cell.convert_xyz_to_integer_abc(out)
        uvw_across = self.cell.convert_xyz_to_integer_abc(across)
        uvw_up     = self.cell.convert_xyz_to_integer_abc(up)
        debugger.print(f"uvw_out   : {out} to {uvw_out}")
        debugger.print(f"uvw_across: {across} to {uvw_across}")
        debugger.print(f"uvw_up    : {up} to {uvw_up}")
        self.current_uvw_orientation = (uvw_up, uvw_across, uvw_out)
        self.current_orientation = self.orientation_definitions[self.orientation]
        self.current_phase = int(self.number_of_phases / 2)
        self.update()

    def define_surface_orientations(self, cell, hkl, uvw):
        """Set the up orientations for the surface hkl.

        hkl defines the surface as (hkl)
        uvw defines a normal to the surface (hkl), [uvw] is made orthogonal to the (hkl) normal
        if (hkl) is (000) [uvw] is used to define the unique direction

        Parameters
        ----------
        cell : unitCell
            The cell used to define the orientations with respect to a, b, c
        hkl : a tuple of 3 integers
            The values of h, k and l
        uvw : a tuple of 3 integers
            The values of u, v and w


        Modifies
        --------
        orientation_definitions

        """        
        debugger.print("define_surface_orientations", cell,hkl,uvw)
        if cell is None: 
            return
        self.cell = cell
        hklsum = abs(hkl[0]) + abs(hkl[1]) + abs(hkl[2])
        if hklsum > 0:
            #
            # Choose the out axis, by converting hkl to xyz
            #
            s = cell.convert_hkl_to_xyz(hkl)
            s = s / np.linalg.norm(s)
            debugger.print("s    :",s)
            #
            # Make sure that the [uvw] direction is orthogonal to (hkl)
            #
            uvw_xyz = cell.convert_abc_to_xyz(uvw)
            uvw_xyz = uvw_xyz / np.linalg.norm(uvw_xyz)
            sdot = np.dot(s, uvw_xyz)
            if abs(sdot) > 0.999:
                #
                # [uvw] and (hkl) are too close to each other
                # Find the lattice vector that is least like the surface normal
                #
                sdot_min = 1.0E12
                i_min = -1
                for i,abc in enumerate(cell.lattice):
                    abc = abc / np.linalg.norm(abc)
                    sdot = np.dot(s,abc)
                    debugger.print(f"abc {abc}, sdot {sdot} ")
                    if abs(sdot) < sdot_min:
                        sdot_min = sdot
                        i_min    = i
                        debugger.print("sdot_min    :",sdot_min)
                        debugger.print("i_min       :",i_min)
                    #
                    # Choose the up axis, by projecting out the surface normal
                    #
                    s_up = cell.lattice[i_min] - np.dot(s,cell.lattice[i_min])*s
            else:
                s_up = uvw_xyz - np.dot(s,uvw_xyz)*s
        else:
            #
            # hkl has been set to (000) so [uvw] becomes the dominant vector
            #
            s = cell.convert_abc_to_xyz(uvw)
            s = s / np.linalg.norm(s)
            #
            # Find the lattice vector that is least like the surface normal
            #
            sdot_min = 1.0E12
            i_min = -1
            for i,abc in enumerate(cell.lattice):
                abc = abc / np.linalg.norm(abc)
                sdot = np.dot(s,abc)
                debugger.print(f"abc {abc}, sdot {sdot} ")
                if abs(sdot) < sdot_min:
                    sdot_min = sdot
                    i_min    = i
                    debugger.print("sdot_min    :",sdot_min)
                    debugger.print("i_min       :",i_min)
            #
            # Choose the up axis, by projecting out the dominant vector
            #
            s_up = cell.lattice[i_min] - np.dot(s,cell.lattice[i_min])*s
        #
        # Choose the across axis by find the normal to s and s_up
        #
        s_across = np.cross(s,s_up)
        s_up     = s_up / np.linalg.norm(s_up)
        s_across = s_across / np.linalg.norm(s_across)
        self.orientation_definitions["s"] = (s_up, s_across, s)
        self.orientation_definitions["S"] = (s, s_across, s_up)
        debugger.print("s       :",s)
        debugger.print("s_up    :",s_up)
        debugger.print("s_across:",s_across)
        return

    def define_orientations(self, cell):
        """Set the up orientations for the cell.

                             Points out    Points up  Points across

                   "x" means     x-axis       y-axis         z-axis
                   "y" means     y-axis       z-axis         x-axis
                   "z" means     z-axis       x-axis         y-axis
                   "a" means     a-axis       axb            ax(axb)
                   "b" means     b-axis       bxc            bx(bxc)
                   "c" means     c-axis       cxa            cx(cxa)
                   "a*" means    a*-axis      a*xb*          a*x(a*xb*)
                   "b*" means    b*-axis      b*xc*          b*x(b*xc*)
                   "c*" means    c*-axis      c*xa*          c*x(c*xa*)

        A capital letter indicates that the up and out designations are reversed


        Parameters
        ----------
        cell : unitCell
            The cell used to define the orientations with respect to a, b, c

        Modifies
        --------
        orientation_definitions

        """        
        debugger.print("define_orientations", cell)
        a = cell.lattice[0] / np.linalg.norm(cell.lattice[0])
        b = cell.lattice[1] / np.linalg.norm(cell.lattice[1])
        c = cell.lattice[2] / np.linalg.norm(cell.lattice[2])
        a_up = np.cross(a,b)
        b_up = np.cross(b,c)
        c_up = np.cross(c,a)
        a_across = np.cross(a,a_up)
        b_across = np.cross(b,b_up)
        c_across = np.cross(c,c_up)
        astar = cell.reciprocal_lattice.T[0] / np.linalg.norm(cell.reciprocal_lattice.T[0])
        bstar = cell.reciprocal_lattice.T[1] / np.linalg.norm(cell.reciprocal_lattice.T[1])
        cstar = cell.reciprocal_lattice.T[2] / np.linalg.norm(cell.reciprocal_lattice.T[2])
        astar_up = np.cross(astar,bstar)
        bstar_up = np.cross(bstar,cstar)
        cstar_up = np.cross(cstar,astar)
        astar_across = np.cross(astar,astar_up)
        bstar_across = np.cross(bstar,bstar_up)
        cstar_across = np.cross(cstar,cstar_up)
        self.orientation_definitions["a"]  = (a_up, a_across, a)
        self.orientation_definitions["a*"] = (astar_up, astar_across, astar)
        self.orientation_definitions["A"]  = (a, a_across, a_up)
        self.orientation_definitions["A*"] = (astar, astar_across, astar_up)
        self.orientation_definitions["b"]  = (b_up, b_across, b)
        self.orientation_definitions["b*"] = (bstar_up, bstar_across, bstar)
        self.orientation_definitions["B"]  = (b, b_across, b_up)
        self.orientation_definitions["B*"] = (bstar, bstar_across, bstar_up)
        self.orientation_definitions["c"]  = (c_up, c_across,c )
        self.orientation_definitions["c*"] = (cstar_up, cstar_across,cstar )
        self.orientation_definitions["C"]  = (c, c_across, c_up)
        self.orientation_definitions["C*"] = (cstar, cstar_across, cstar_up)
        debugger.print("a       :",a)
        debugger.print("a_up    :",a_up)
        debugger.print("a_across:",a_across)
        debugger.print("b       :",b)
        debugger.print("b_up    :",b_up)
        debugger.print("b_across:",b_across)
        debugger.print("c       :",c)
        debugger.print("c_up    :",c_up)
        debugger.print("c_across:",c_across)
        return

    def save_movie(self, filename):
        """Save the current state as a movie file.

        This method captures snapshots of the current state at different phases,
        processes them, and compiles them into a movie file.

        Parameters
        ----------
        filename : str
            The filename (including path) where the movie should be saved.

        Returns
        -------
        None

        Notes
        -----
        Relies on `imageio` for writing the movie file and handling images.
        If a timer is present, it's temporarily stopped to ensure snapshots are
        taken consistently, and then restarted once the movie is saved.
        Temporary storage is used for intermediate images

        """        
        import imageio
        # imageio.plugins.ffmpeg.download()
        debugger.print("save_movie", filename)
        if self.timer is not None:
            self.timer.stop()
        writer  = imageio.get_writer(filename, mode="I", fps=24)
        tmpdir  = os.path.dirname(filename)
        tmpfile = os.path.join(tmpdir,".snapshot.png")
        debugger.print("save_movie", filename)
        for _i in range(0,2*self.number_of_phases):
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
            writer.append_data(image)
        writer.close()
        os.remove(tmpfile)
        if self.timer is not None:
            self.timer.start()

    def snapshot(self,filename):
        """Capture and save a snapshot of the current framebuffer to a file.

        Parameters
        ----------
        filename : str
            The name of the file where the snapshot will be saved.

        Returns
        -------
        None

        Examples
        --------
        Assuming an instance `display` of the class containing this method:

        ```python
        display.snapshot('screenshot.png')
        ```

        This will grab the current framebuffer content and save it as 'screenshot.png'.

        """        
        debugger.print("snapshot", filename)
        image = self.grabFramebuffer()
        image.save(filename)

    def translate(self,x,y):
        """Translate an object by modifying its current position.

        Parameters
        ----------
        x : float
            The distance to translate along the x-axis.
        y : float
            The distance to translate along the y-axis.

        Returns
        -------
        None

        Notes
        -----
        This method updates the current position of the object by adding the specified distances along the x and y axes, respectively. It does not return a value.

        This function involves a call to `glTranslatef`, which applies a translation transformation to the current matrix. Ensure that the appropriate OpenGL context is made current before calling this function by using `self.myMakeCurrent()`.

        """        
        debugger.print("translate ",x,y)
        self.myMakeCurrent()
        glTranslatef(x, y, 0.0)

    def wheelEvent(self, event):
        """Handle a wheel event to zoom in or out.

        Parameters
        ----------
        event : QWheelEvent
            The wheel event containing information about the wheel movement.

        Returns
        -------
        None

        Notes
        -----
        This function first logs the occurrence of a wheel event. It then sets
        the current context to this object before determining the direction and
        intensity of the zoom based on the vertical angle delta of the wheel event.
        Lastly, it refreshes the display to reflect any changes.

        """        
        debugger.print("Wheel event " )
        self.myMakeCurrent()
        zoom = event.angleDelta().y()
        self.zoom(zoom)
        self.update()

    def mousePressEvent(self, event):
        """Handle the mouse press event.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event instance containing data of the mouse press.

        Returns
        -------
        None

        """        
        self.xAtPress = event.x()
        self.yAtPress = event.y()
        #print('mouse press',self.xAtPress, self.yAtPress)

    def mouseReleaseEvent(self, event):
        """Handle the mouse release event.

        Parameters
        ----------
        event : QMouseEvent
            The event object containing information about the mouse release event.

        Returns
        -------
        None

        Notes
        -----
        This method records the x and y coordinates of the mouse at the moment it is released.

        """        
        self.xAtRelease = event.x()
        self.yAtRelease = event.y()
        #print('mouse release',self.xAtRelease, self.yAtRelease)

    def zoom(self, zoom):
        """Zoom in or out on an object.

        Parameters
        ----------
        zoom : int or float
            The zoom level. Positive values zoom in, negative values zoom out.

        Returns
        -------
        None

        Notes
        -----
        - This function requires a current context and will make itself the current one before applying transformations.
        - Updates the object's state after zooming.

        Uses `glScalef` from OpenGL for scaling transformations and assumes that `debugger` and `glScalef` are accessible in the current context, with `debugger.print` used for debug logging.

        """        
        debugger.print("zoom ", zoom)
        self.myMakeCurrent()
        zoom_factor = 1.06 if zoom > 0 else 0.94
        debugger.print("zoom factor", zoom_factor)
        glScalef(zoom_factor, zoom_factor, zoom_factor)
        self.update()

    def mouseMoveEvent(self, event):
        """Handle the mouse move events to perform various actions.

        Parameters
        ----------
        event : QMouseEvent
            The mouse event that contains information about the state of the mouse.

        Returns
        -------
        None

        Notes
        -----
        This method checks if the left or middle mouse button is pressed and applies different transformations (zoom, rotate, or translate) based on the mouse movement and whether certain modifier keys (Shift or Control) are pressed.
        - If the left mouse button is pressed:
        - With no modifiers or with Shift or Control modifier pressed: Performs rotation.
        - With Shift or Control modifier pressed: Performs zoom.
        - If the middle mouse button is pressed: Performs translation.

        Mouse coordinates are updated at the end of the method to reflect the current position.

        """        
        self.xAtMove = event.x()
        self.yAtMove = event.y()
        buttons = event.buttons()
        modifiers = event.modifiers()
        if buttons & Qt.LeftButton:
            debugger.print("Mouse event - left button")
            if modifiers & Qt.ShiftModifier or modifiers & Qt.ControlModifier:
                # handle zoom
                xzoom = +1.0*(self.xAtMove - self.xAtPress)
                yzoom = -1.0*(self.yAtMove - self.yAtPress)
                self.zoom(xzoom+yzoom)
            else:
                # handle ordinary rotation
                xrotate =-(self.xAtMove - self.xAtPress)
                yrotate = (self.yAtMove - self.yAtPress)
                self.moleculeRotate(0.3,yrotate,xrotate,0)
        elif buttons & Qt.MidButton:
            debugger.print("Mouse event - mid button")
            xshift = -0.02 * (self.xAtMove - self.xAtPress)
            yshift = -0.02 * (self.yAtMove - self.yAtPress)
            up,across,out = self.current_orientation
            shifted = xshift*np.array(across)+yshift*np.array(up)
            self.translate(shifted[0], shifted[1])
            self.update()
        debugger.print("Mouse event - xy", self.xAtPress,self.yAtPress)
        self.xAtPress = self.xAtMove
        self.yAtPress = self.yAtMove

    def stopAnimation(self):
        """Stop the ongoing animation.

        This method stops the animation by stopping the timer associated with it.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This function assumes that `self.timer` is an instance with a `stop` method, which is called to halt any ongoing animation. If `self.timer` is `None`, indicating no animation is currently active, this method does nothing.

        """        
        debugger.print("stopAnimation")
        if self.timer is not None:
            self.timer.stop()

    def startAnimation(self):
        """Start or restarts an animation timer.

        This method checks if an animation timer is already running.
        If so, it stops the current timer. If not, it initializes and starts a new timer
        with the interval specified in `self.timer_interval`. When the timer times out,
        it triggers `self.timeoutHandler`.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Notes
        -----
        - `self.timer` is expected to be an attribute of the instance which either holds
          a QTimer object or None.
        - This method directly modifies `self.timer` to either stop the current timer
          or create, set up, and start a new timer.
        - This method is designed for use within a class that manages animations or
          timed events, utilizing a QTimer from the PyQt or PySide frameworks.

        """        
        debugger.print("startAnimation")
        if self.timer is not None:
            self.timer.stop()
        else:
            self.timer = QTimer(self)
            self.timer.timeout.connect(self.timeoutHandler)
            self.timer.setInterval(self.timer_interval)
        self.timer.start()

    def paintGL(self):
        """Render OpenGL graphics for the current scene.

        This method prepares the graphics context, clears the buffer, sets the background color,
        and draws the current scene which may include spheres, cylinders, texts and arrows based on the object's state.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Notes
        -----
        Requires access to OpenGL context and assumes that necessary OpenGL libraries
        and the debugger object are available in the scope. `paintGL` is typically called
        by the OpenGL rendering loop.

        Attributes Modified
        -------------------
        - Utilizes several attributes of the instance (`self`), including `background_colour`,
          `rotation_centre`, and `show_arrows`, to determine how to render the scene.
        - The OpenGL context is modified by setting the clear color and transforming the modelview matrix.

        """        
        debugger.print("paintGL")
        glMatrixMode(GL_MODELVIEW)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        self.background_colour = np.array(self.viewerTab.settings["Background colour"])/255.0
        glClearColor(*self.background_colour)
        glPushMatrix()
        glTranslatef(-self.rotation_centre[0],-self.rotation_centre[1],-self.rotation_centre[2] )
        self.drawSpheres()
        self.drawCylinders()
        self.drawTexts()
        if self.show_orientation:
            self.drawUVWInfo()
        if self.show_arrows:
            self.drawArrows()
        glPopMatrix()

    def drawUVWInfo(self):
        """Draw the uvw info on the screen.

        The uvw information (up, across, out of the screen)

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("drawUVWInfo")
        uvw_up, uvw_across, uvw_out = self.current_uvw_orientation
        uvw_string = f"Normal to screen = {uvw_out},  screen vertical = {uvw_up}"
        self.renderText(10,10,0,uvw_string,screen_coordinates=True)
        debugger.print("drawUVWInfo", uvw_string)
        return

    def drawSpheres(self):
        """Draw spheres based on the object's current phase and spheres attributes.

        This method iterates through the spheres assigned to the current phase of the invoking object. For each sphere, it sets the material properties based on the sphere's color and predefined material factors. Then, it translates the drawing location to the sphere's position and draws the sphere using OpenGL functions.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Notes
        -----
        Requires OpenGL (PyOpenGL) for rendering. The method modifies the OpenGL state to draw spheres at specified locations with given radii and colors. Assumes the presence of attributes like `spheres`, `current_phase`, `diffuseMaterialFactor`, `ambientMaterialFactor`, `specularLightFactor`, and methods or attributes for rendering parameters such as `quadric`, `sphere_slices`, `sphere_stacks`.

        """        
        debugger.print("drawSpheres")
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

    def drawTexts(self):
        """Draw the texts stored in the object.

        This method iterates through the list of text objects for the current phase, 
        computes the required material colors based on predefined factors, and renders 
        each text using OpenGL commands.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("drawTexts")
        if len(self.texts) == 0:
            return
        for text,_color,_size,pos in self.texts[self.current_phase]:
            self.renderText(pos[0],pos[1],pos[2],text)

    def renderText(self, x, y, z, string, screen_coordinates=False):
        """Draw text in the widget.

        Draw the text in string at the position starting x, y, z

        Parameters
        ----------
        x, y, z : floats
            The position at which to start drawning
        string : str
            The text
        screen_coordinates : bool
            Optional parameter by default the coordinates are model coordinates
            If true the screen coordinates are used

        Returns
        -------
        None

        """        
        height = self.height()

        model = glGetDoublev(GL_MODELVIEW_MATRIX)
        proj  = glGetDoublev(GL_PROJECTION_MATRIX)
        view  = glGetIntegerv(GL_VIEWPORT)
        if not screen_coordinates:
            textPosX, textPosY, textPosZ = gluProject(x, y, z, model, proj, view)
            textPosY = height - textPosY
        else:
            textPosX = x
            textPosY = height-y
        glPushAttrib(GL_ALL_ATTRIB_BITS)
        painter = QPainter()
        painter.begin(self)
        painter.setPen(Qt.yellow)
        painter.setFont(QFont("Helvetica", 14))
        painter.setRenderHints(QPainter.Antialiasing | QPainter.TextAntialiasing)
        painter.drawText(int(textPosX), int(textPosY), string)
        painter.end()
        glPopAttrib()
        return


    def drawCylinders(self):
        """Draw the cylinders stored in the object.

        This method iterates through the list of cylinder objects for the current phase, 
        computes the required material colors based on predefined factors, and renders 
        each cylinder using OpenGL commands.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Notes
        -----
        Each cylinder is represented by an object with attributes such as colour, radius,
        position1, position2, height, angle, and rotation. These attributes are used to 
        determine the cylinder's appearance and placement in the 3D space.

        This function assumes that OpenGL context is correctly set up, and global variables 
        such as `self.diffuseMaterialFactor`, `self.ambientMaterialFactor`, `self.specularLightFactor`,
        `self.white`, `self.glintMaterialFactor`, `self.quadric`, `self.cylinder_slices`, and 
        `self.cylinder_stacks` are predefined and configured.

        It is also assumed that `self.cylinders` is a dictionary where each key corresponds to a 
        different phase, and the value is a list of cylinder objects for that phase. The current 
        phase is determined by `self.current_phase`.

        """        
        debugger.print("drawCylinders")
        if len(self.cylinders) == 0:
            return
        for cylinder in self.cylinders[self.current_phase]:
            col    = cylinder.colour
            rad    = cylinder.radius
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
        """Draw the designated arrows with specified properties.

        This method draws a series of arrows based on the properties stored within instances of the class. The arrows are displayed on the viewer associated with `self`. Each arrow has an associated sphere from which it originates. These spheres have positional data that dictate where the arrow is drawn. Furthermore, arrows have several visual properties such as color and size that are applied when rendering. Error handling for absent arrows is included, halting the function early if no arrows are to be drawn.

        Parameters
        ----------
        None

        Returns
        -------
        None


        Notes
        -----
        - The arrow drawing process involves translating and rotating the OpenGL context to the appropriate position and orientation before drawing cylinders to represent the arrows.
        - Material properties such as diffuse, ambient, and specular colors are calculated from the arrow's color and applied to the OpenGL material.
        - Arrows consist of a cylindrical shaft and a conical tip, whose dimensions and orientations are based on the properties of each arrow (e.g., direction, length, angle, rotation, colour, and radius).
        - This method depends on the `self.arrows` list and the corresponding `self.spheres` for the current phase stored in `self.current_phase`. It leverages OpenGL functions and assumes that `self.viewerTab.settings['Arrow colour']`, along with other properties like `self.diffuseMaterialFactor`, `self.ambientMaterialFactor`, `self.specularLightFactor`, and `self.glintMaterialFactor` are predefined and correctly configured.

        """        
        debugger.print("drawArrows")
        if len(self.arrows) == 0:
            return
        self.arrow_colour = np.array(self.viewerTab.settings["Arrow colour"])/255.0
        for arrow,sphere in zip(self.arrows,self.spheres[self.current_phase]):
            pos       = sphere.position
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
        """Resizes the OpenGL widget to a new width and height.

        Notes
        -----
        This is typically overridden in a subclass and used to adjust the OpenGL viewport or projection matrix when the widget is resized.

        Parameters
        ----------
        w : int
            The new width of the OpenGL widget.
        h : int
            The new height of the OpenGL widget.

        Returns
        -------
        None

        """        
        debugger.print("resizeGL",w,h)
        self.my_width = w
        self.my_height =h
        # set projection matrix
        self.setProjectionMatrix()

    def initializeGL(self):
        """Initialise OpenGL rendering context.

        This function sets up various OpenGL attributes and states including lighting, depth test, smooth shading, and anti-aliasing features to prepare for rendering in a GL context. It configures the viewport settings based on the associated viewer tab's background color settings.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("initializeGL")
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
        self.background_colour = np.array(self.viewerTab.settings["Background colour"])/255.0
        glClearColor(*self.background_colour)

    def setImageSize(self):
        """Set the maximum image size based on the furthest object from the rotation center.

        This method evaluates the position of spheres and cylinders in the current phase, calculates their distance from a fixed point (rotation center), and updates the `image_size` attribute to the maximum distance found. It ensures the all objects will fit within the rendered image. Finally, it updates the projection matrix and logs the new image size.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        The `self` parameter represents the instance of the class in which the function is called. This function relies on the `math.sqrt`, `np.dot` for calculations, and expects `self.spheres`, `self.cylinders`, `self.current_phase`, `self.rotation_centre`, and `self.setProjectionMatrix` to be defined.

        """        
        maxsize = 0.0
        for sphere in self.spheres[self.current_phase]:
            pos = sphere.position
            vec = pos - self.rotation_centre
            dist = math.sqrt(np.dot(vec,vec))
            maxsize = max(maxsize, dist)
        if maxsize > 0:
            self.image_size = maxsize
        else:
            self.image_size = 10.0
        self.setProjectionMatrix()
        debugger.print("setImageSize",self.image_size)

    def setProjectionMatrix(self):
        """Set the projection matrix for rendering.

        This method configures the projection matrix based on the current image size, width, and height attributes of the object. It's intended to set up an orthographic projection suited to the object's dimensions. If any of the required attributes (`image_size`, `my_width`, or `my_height`) are not set, the method will exit early. The method also sets the view to the middle phase of an unspecified number of total phases and resets the model view matrix.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method makes several OpenGL calls to configure the projection and model view matrices and assumes that `self.myMakeCurrent()` makes the required OpenGL context current.
        It also updates the rotation with `self.matrix` and sets `self.current_phase` based on the total number of phases.

        """        
        if self.image_size is None or self.my_width is None or self.my_height is None:
            return
        self.myMakeCurrent()
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        orthox = 1.1 * self.image_size * self.my_width  / min(self.my_width,self.my_height)
        orthoy = 1.1 * self.image_size * self.my_height / min(self.my_width,self.my_height)
        orthoz = 1.1 * self.image_size * max(self.my_width,self.my_height)
        debugger.print("projection", orthox, orthoy, orthoz)
        glOrtho(-orthox, orthox, -orthoy, orthoy, -orthoz, orthoz)
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        glMultMatrixf(self.matrix)
        # reset the current phase to the centre of the phases
        self.current_phase = int(self.number_of_phases / 2)
        debugger.print("set projection matrix ortho", orthox, orthoy, orthoz)
        debugger.print("set projection matrix image_size", self.image_size)
        debugger.print("set projection matrix current_phase", self.current_phase)

    def defineLights(self):
        """Define light configurations for the viewer.

        This method sets up the lighting for a 3D viewer by configuring ambient, diffuse, and specular light properties, as well as light positions. It also controls which lights are enabled or disabled based on the object's `light_switches` attribute.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        - `self.light_switches` needs to be initially set or inherited from `self.viewerTab.light_switches` if `self.light_switches` is None.
        - Lighting is enabled globally with `glEnable(GL_LIGHTING)`.
        - For each light defined by the `self.light_positions`, `self.light_switches`, and `self.lights` attributes, this method sets up:
          - Ambient light properties with `glLightfv(light, GL_AMBIENT, [...])`.
          - Diffuse light properties with `glLightfv(light, GL_DIFFUSE, [...])`.
          - Specular light properties with `glLightfv(light, GL_SPECULAR, [...])`.
          - The position of each light with `glLightfv(light, GL_POSITION, position)`.
        - If `self.linearAttenuation` is True, linear attenuation is set for each light with `glLight(light, GL_LINEAR_ATTENUATION, 1.0)`.
        - Each light in `self.lights` is either enabled or disabled based on its corresponding status in `self.light_switches`.

        """        
        debugger.print("Define Lights")
        self.myMakeCurrent()
        if self.light_switches is None:
            self.light_switches = self.viewerTab.light_switches
            glEnable(GL_LIGHTING)
            for position,light in zip(self.light_positions,self.lights):
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
        """Set the rotation centre to a specified position.

        Parameters
        ----------
        pos : array_like
            The position to set the rotation centre. Expected to be an iterable
            with numerical values that can be converted to a NumPy array.

        Returns
        -------
        None

        """        
        debugger.print("set rotation centre",pos)
        self.rotation_centre = np.array(pos)

    def createArrays(self, nphases):
        """Create or reinitialize arrays for spheres and cylinders for each phase.

        This method creates or clears and then repopulates `self.spheres` and `self.cylinders` with an empty `deque` for each phase. It also sets `self.number_of_phases` to the given `nphases` and initializes `self.current_phase` to the midpoint of the range of phases (rounded down if `nphases` is an odd number).

        Parameters
        ----------
        nphases : int
            The number of phases for which to create arrays.

        Returns
        -------
        None

        Notes
        -----
        - Assumes `self.spheres`, `self.texts` and `self.cylinders` are already defined as lists.
        - It uses the `clear` method, which is available for lists in Python 3.3 and later.

        """        
        debugger.print("createArrays")
        #  Create an empty list for each phase
        self.spheres.clear()
        self.cylinders.clear()
        self.texts.clear()
        for _i in range(nphases):
            self.spheres.append(deque())
            self.cylinders.append(deque())
            self.texts.append(deque())
        self.number_of_phases = nphases
        self.current_phase = int(self.number_of_phases / 2)

    def deleteSpheres(self):
        """Delete all spheres from the current context.

        This method clears the collection of spheres.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("deleteSpheres")
        self.spheres.clear()

    def deleteTexts(self):
        """Delete all texts from the current context.

        This method clears the collection of texts.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("deleteTexts")
        self.texts.clear()

    def deleteCylinders(self):
        """Delete all cylinder objects from the storage.

        This method clears the list that stores cylinder objects, effectively
        removing all the cylinders managed by the instance.

        Parameters
        ----------
        None

        Returns
        -------
        None

        Notes
        -----
        This method should be used with caution as it irreversibly deletes all 
        cylinder data from the instance.

        """        
        debugger.print("deleteCylinders")
        self.cylinders.clear()

    def deleteArrows(self):
        """Delete all arrows from a collection.

        This method clears the collection of arrows by removing all elements from it.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """        
        debugger.print("deleteArrows")
        #self.arrows = []
        self.arrows.clear()

    def addArrows(self, colour, radius, direction, length, phase=0):
        # There is no phase requirement for arrows - they are just displacements
        """Add an arrow to the list of arrows.

        Parameters
        ----------
        colour : str
            The color of the arrow.
        radius : int or float
            The radius from the center where the arrow will be placed.
        direction : int or float
            The direction in which the arrow points.
        length : int or float
            The length of the arrow.
        phase : int or float, optional
            The phase of the arrow, defaults to 0.

        Returns
        -------
        None

        """        
        self.arrows.append( Arrow(colour, radius, direction, length) )

    def addText(self, string, colour, size, pos, phase=0):
        """Add a text object to the specified phase collection of texts.

        Parameters
        ----------
        string : str
            The content of the text.
        colour : str
            The color of the text.
        size : float
            The size of the text
        pos : tuple
            The starting position of the text (x, y, z coordinates).
        phase : int, optional
            The phase to which the text should be added. Default is 0.

        Returns
        -------
        None

        """        
        self.texts[phase].append( (string,colour,size,pos) )

    def addCylinder(self, colour, radius, pos1, pos2, phase=0):
        """Add a Cylinder object to the specified phase collection of cylinders.

        Parameters
        ----------
        colour : str
            The color of the cylinder.
        radius : float
            The radius of the cylinder.
        pos1 : tuple
            The starting position of the cylinder (x, y, z coordinates).
        pos2 : tuple
            The ending position of the cylinder (x, y, z coordinates).
        phase : int, optional
            The phase to which the cylinder should be added. Default is 0.

        Returns
        -------
        None

        """        
        self.cylinders[phase].append( Cylinder(colour,radius,pos1,pos2) )

    def addSphere(self, colour, radius, pos, phase=0):
        """Add a sphere to a collection within an object.

        Parameters
        ----------
        colour : Any
            The colour of the sphere, can vary depending on representation (e.g., string, tuple of RGB values).
        radius : float
            The radius of the sphere.
        pos : tuple
            Position of the sphere, usually a tuple of coordinates (x, y, z).
        phase : int, optional
            The specific phase (or collection) to which the sphere should be added. Defaults to 0.

        Returns
        -------
        None

        """        
        self.spheres[phase].append( Sphere(colour,radius,pos) )

class Sphere:
    """A class representing a sphere.

    Parameters
    ----------
    colour : list of int
        The RGB colour of the sphere as a list of three integers ranging from 0 to 255.
    radius : float
        The radius of the sphere.
    position : list of float
        The position of the sphere in 3D space as a list of three coordinates (x, y, z).

    Attributes
    ----------
    colour : ndarray
        The normalized RGB colour of the sphere, with each element ranging from 0.0 to 1.0.
    radius : float
        The radius of the sphere.
    position : ndarray
        The position of the sphere in 3D space as a NumPy array with three elements (x, y, z).

    """

    def __init__(self, colour, radius, position):
        """Initialize a new object with colour, radius, and position.

        Parameters
        ----------
        colour : list or np.ndarray
            A list or an array of RGB values (0 to 255) for the colour.
        radius : float
            The radius of the object.
        position : list or np.ndarray
            A list or an array representing the position in 2D or 3D space.

        Notes
        -----
        - The `colour` parameter will be normalized to a range between 0 and 1 internally.
        - The `position` parameter is expected to match the dimensionality of the system (i.e., 2D or 3D).

        """        
        self.colour= np.array(colour)/255.0
        self.radius = radius
        self.position = np.array(position)

class Cylinder:
    """A class representing a 3D cylinder.

    Parameters
    ----------
    colour : array_like
        The colour of the cylinder, expected to be a list or array of three elements (R, G, B) where each value is in the range 0-255.
    radius : float
        The radius of the cylinder.
    pos1 : array_like
        The 3D coordinates (as a list or array) of one end of the cylinder's central axis.
    pos2 : array_like
        The 3D coordinates (as a list or array) of the other end of the cylinder's central axis.

    Attributes
    ----------
    colour : ndarray
        The normalized colour of the cylinder (R, G, B) with values ranging from 0.0 to 1.0.
    radius : float
        The radius of the cylinder.
    position1 : ndarray
        The 3D position of one end of the cylinder's central axis.
    position2 : ndarray
        The 3D position of the other end of the cylinder's central axis.
    angle : float
        The angle in degrees between the cylinder's central axis and the z-axis.
    height : float
        The height of the cylinder, determined as the distance between `pos1` and `pos2`.
    rotation : ndarray
        The 3D vector representing the axis of rotation to align the cylinder with the z-axis (used together with `angle` for setting the cylinder's orientation in 3D space).

    """

    def __init__(self, colour, radius, pos1, pos2):
        """Initialize a new instance with specified color, radius, and positions.

        Parameters
        ----------
        colour : array_like
            RGB values of the color, each component must be between 0 and 255.
        radius : float
            The radius of the instance.
        pos1 : array_like
            The first position vector, described by three coordinates. 
        pos2 : array_like
            The second position vector, also described by three coordinates.

        Notes
        -----
        This method computes a vector from `pos2` to `pos1`, calculates its height,
        and determines the angle of rotation from the vertical axis. The color is
        normalized to be in the range [0, 1]. Positions and rotation are stored as
        NumPy arrays, and properties such as `height`, `angle`, and `rotation` direction
        are calculated for use in positioning and rendering the instance.

        """        
        if pos1[2] > pos2[2]:
            a = np.array(pos1)
            b = np.array(pos2)
        else:
            b = np.array(pos1)
            a = np.array(pos2)
        self.colour= np.array(colour)/255.0
        self.radius = radius
        self.position1 = a
        self.position2 = b
        z = np.array( [0.0, 0.0, 1.0] )
        p = a - b
        height = math.sqrt(np.dot(p,p))
        p = p / height
        t = np.cross(z, p)
        angle = 180.0 / np.pi * math.acos( np.dot(z,p) )
        self.angle = angle
        self.height = height
        self.rotation = t

class Arrow:
    """A class to represent a 3D arrow for graphical display.

    This class encapsulates properties of a 3D arrow including its color,
    size (radius), direction, and scaling factor. It uses numpy for
    calculations to set up the appropriate direction, height, and rotation
    angle for the arrow based on the provided parameters.

    Parameters
    ----------
    colour : list of int
        The RGB values of the arrow colour, where each value is in the range [0, 255].
    radius : float
        The radius of the arrow's base.
    direction : list or numpy.ndarray
        A 3-element list or array defining the direction of the arrow in 3D space.
    scale : float
        A scaling factor to adjust the size of the arrow.

    Attributes
    ----------
    colour : numpy.ndarray
        The normalized RGB colour of the arrow, with each channel value in the range [0.0, 1.0].
    angle : float
        The angle between the arrow direction and the z-axis, in degrees.
    height : float
        The length of the arrow scaled by the `scale` parameter.
    rotation : numpy.ndarray
        A 3-element array representing the rotation axis for aligning the arrow with its direction.
    direction : numpy.ndarray
        A normalized 3-element array representing the direction of the arrow.
    radius : float
        The radius of the arrow's base as provided during instantiation.

    """

    def __init__(self, colour, radius, direction, scale):
        """Initialize a custom object with given parameters.

        Parameters
        ----------
        colour : list or np.array
            A 3-element list or array representing the RGB values of the colour. Each element in the range 0-255.
        radius : float
            The radius of the object.
        direction : list or np.array
            A 3-element list or array indicating the direction of the object in 3D space.
        scale : float
            A scaling factor for the object's height.

        Notes
        -----
        This initialization function converts the `colour` from an integer RGB format (0-255 per channel) to 
        a normalized format (0-1.0 per channel). It calculates the rotation and direction vectors based on 
        the input `direction` and `scale`, prepares the object for further use.

        This function uses NumPy operations for calculations. It assumes the direction 
        is not parallel to the z-axis to avoid division by zero in the cross product.

        The angle calculation is done in degrees and considers the case when the direction vector 
        is almost aligned with the z-axis by adding `1e-08` to the height to avoid division by zero.

        Attributes
        ----------
        colour : np.array
            The normalized RGB colour of the object.
        radius : float
            The radius of the object.
        direction : np.array
            The normalized direction vector of the object.
        scale : float
            The scale factor for the object's height.
        angle : float
            The angle between the direction vector and the z-axis, in degrees.
        height : float
            The scaled height of the object based on the input direction and scale.
        rotation : np.array
            The rotation vector of the object.

        """        
        # Adds a displacement to each atom
        self.colour= np.array(colour)/255.0
        p = np.array(direction)
        z = np.array( [0.0, 0.0, 1.0] )
        t = np.cross(z, p)
        height = math.sqrt(np.dot(p,p)) + 1.0E-8
        angle = 180.0 / np.pi * math.acos( np.dot(z,p) / height )
        p = p / height
        self.angle = angle
        self.height = height*scale
        self.rotation = t
        self.direction = p
        self.radius = radius

