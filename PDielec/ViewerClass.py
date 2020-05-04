from traits.api import HasTraits, Button, Range, Float, Int, Instance, on_trait_change, Property
from traitsui.api import View, Item, Group, HGroup, RangeEditor, Label, spring
from mayavi.core.api import PipelineBase
from mayavi.core.ui.api import MayaviScene, SceneEditor, MlabSceneModel
import numpy as np
from mayavi import mlab
from PDielec.Constants import covalent_radii, elemental_colours

class ViewerClass(HasTraits):
    number_of_modesm1 = 0
    maxdisplacement = 0.3
    gui_mode_selection    = Int(label="Mode")
    gui_maxdisplacement   = Float(maxdisplacement,label="Displacement")
    gui_savepicture_button = Button('Snapshot')
    gui_animate_button    = Button('Animate')
    gui_exit_button       = Button('Exit')
    # gui_exit_button.label = '?????'
    #from objbrowser import browse
    #browse(gui_exit_button)
    selected_mode_frequency = Property(Float, depends_on = ['gui_mode_selection'], label="Fequency")
#
    scene = Instance(MlabSceneModel, ())
    plot = Instance(PipelineBase)
    view = View(Item('scene', editor=SceneEditor(scene_class=MayaviScene),
                     height=250, width=300, show_label=False,resizable=True),
                Group(
                    Item('gui_mode_selection',
                         editor =  RangeEditor(low=0, high=number_of_modesm1, high_name='number_of_modesm1',format='%d',mode='auto')
                        ),
#jk                    Label('Select mode number for display of displacements'),
                    Item('selected_mode_frequency',style='readonly'),
                    'gui_maxdisplacement',
                    label='Mode selector',
                    show_labels = True,

                ),
                HGroup(
#                    'gui_savepicture_button',
                    Item('_'),
                    spring,
                    'gui_animate_button',
                    spring,
                    'gui_exit_button',
                    spring,
                    show_labels = False,
                ),
		resizable=True,
           )


    @on_trait_change('scene.activated')
    def update_plot(self):
        #jk print('update_plot',self.gui_mode_selection)
        if self.plot is None:
            #jk print('plot is None')
            self.draw()
        else:
            #jk print('plot None')
            pass

    @on_trait_change('gui_animate_button')
    def update_plot_on_animate_button(self):
        #print('Gui animate button has been pressed ')
        if self.animation:
           self.animation = False
           self.draw()
        else:
           self.animation = True
           self.animate()

    @on_trait_change('gui_savepicture_button')
    def update_plot_on_savepicture(self):
        print('Gui savepicture button has been pressed ')
        filename = 'snapshot_'+str(self.snapshot_number)+'.pdf'
        print(filename)
        mlab.savefig(filename, figure=mlab.gcf(),magnification=3)


    @on_trait_change('gui_maxdisplacement')
    def update_plot_on_maxdisplacement(self):
        #jk print('Gui maxdisplacement button has been pressed ')
        self.maxdisplacement = self.gui_maxdisplacement
        self.calculatePhasePositions()

    @on_trait_change('gui_exit_button')
    def update_plot_on_exit_button(self):
        #print('Gui exit button has been pressed ')
        exit()

    @on_trait_change('gui_mode_selection')
    def update_plot_on_mode_change(self):
        #print('mode has changed ',self.gui_mode_selection)
        self.mode = self.gui_mode_selection
        self.calculatePhasePositions()
        self.refreshDisplacements()

    def __init__(self, unit_cell, normal_modes, frequencies):
        super(ViewerClass, self).__init__()
        self.quiverBonds = False
        self.snapshot_number = 0
        self.unit_cell = unit_cell
        self.normal_modes = normal_modes
        self.frequencies = frequencies
        self.bonds = unit_cell.bonds
        self.nbonds = len(self.bonds)
        self.scale_vibrations = 20.0
        self.animation = False
        self.X = self.unit_cell.xyz_coordinates[:,0]
        self.Y = self.unit_cell.xyz_coordinates[:,1]
        self.Z = self.unit_cell.xyz_coordinates[:,2]
        self.natoms = len(self.X)
        self.number_of_modes = len(normal_modes)
        self.number_of_modesm1 = len(normal_modes)-1
        # reorder the displacement info in the normal modes into U,V and W lists
        self.U, self.V, self.W = [], [], []
        for mode,displacements in enumerate(self.normal_modes):
            U = np.zeros_like( self.X )
            V = np.zeros_like( self.X )
            W = np.zeros_like( self.X )
            index = -1
            for i in range(0,len(displacements),3):
                index += 1
                U[index] = displacements[i + 0]
                V[index] = displacements[i + 1]
                W[index] = displacements[i + 2]
            self.U.append(U)
            self.V.append(V)
            self.W.append(W)
        # get the cell edges for the bounding box
        corners,self.cell_edges = self.unit_cell.getBoundingBox()
        self.element_names = self.unit_cell.element_names
        self.colours = []
        self.radii = []
        # set the colours (the colour is a triplet of floats from 0.0 to 1.0)
        for el in self.element_names:
           colour_list = [ colour/255.0 for colour in elemental_colours[el] ]
           colour_tuple = tuple(colour_list)
           self.colours.append( colour_tuple )
           self.radii.append(1.3*covalent_radii[el])
        # We need these glyphs to update the atom and bond positions for animations
        self.glyph_atoms = []
        self.glyph_bonds = []
        self.mode = 0
        self.number_of_phase_steps = 21
        self.phase_direction = 1
        self.phase_index = 0
        self.calculatePhasePositions()
        return

    def calculatePhasePositions(self):
        # we need the number of phase steps to be odd
        if self.number_of_phase_steps%2 == 0:
            self.numbers_of_phase_steps += 1
        U = self.U[self.mode]
        V = self.V[self.mode]
        W = self.W[self.mode]
        maxR = np.amax(U)
        maxR = max(maxR, np.amax(V) )
        maxR = max(maxR, np.amax(W))
        self.scale_vibrations = self.maxdisplacement / maxR
        #print('MaxR is',maxR)
        self.newX    = np.zeros( (self.number_of_phase_steps,self.natoms) )
        self.newY    = np.zeros( (self.number_of_phase_steps,self.natoms) )
        self.newZ    = np.zeros( (self.number_of_phase_steps,self.natoms) )
        self.startX  = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.startY  = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.startZ  = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.dirX    = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.dirY    = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.dirZ    = np.zeros( (self.number_of_phase_steps,self.nbonds) )
        self.bondx   = np.zeros( (self.number_of_phase_steps,self.natoms,2) )
        self.bondy   = np.zeros( (self.number_of_phase_steps,self.natoms,2) )
        self.bondz   = np.zeros( (self.number_of_phase_steps,self.natoms,2) )
        # work out the scaling factor to get the require maximumdisplacement
        sign = +1.0
        phases = []
        small = 1.0E-8
        n2 = int(self.number_of_phase_steps/2)
        delta = 1.0 / n2
        phase = -delta
        for phase_index in range(self.number_of_phase_steps):
            phase += delta*sign
            self.newX[phase_index] = self.X+U*phase*self.scale_vibrations
            self.newY[phase_index] = self.Y+V*phase*self.scale_vibrations
            self.newZ[phase_index] = self.Z+W*phase*self.scale_vibrations
            if phase >= 1.0-small:
                sign = -1.0
                phase += delta*sign
            elif phase <= -1.0+small:
                sign = +1.0
                phase += delta*sign
            phases.append(phase)
        # end for phase_index
        for phase_index,phase in enumerate(phases):
            for ibond,bond in enumerate(self.bonds):
                i,j = bond
                self.bondx[phase_index,ibond,0] = self.newX[phase_index,i]
                self.bondy[phase_index,ibond,0] = self.newY[phase_index,i]
                self.bondz[phase_index,ibond,0] = self.newZ[phase_index,i]
                self.bondx[phase_index,ibond,1] = self.newX[phase_index,j]
                self.bondy[phase_index,ibond,1] = self.newY[phase_index,j]
                self.bondz[phase_index,ibond,1] = self.newZ[phase_index,j]
                self.startX[phase_index,ibond] = self.newX[phase_index,i]
                self.startY[phase_index,ibond] = self.newY[phase_index,i]
                self.startZ[phase_index,ibond] = self.newZ[phase_index,i]
                self.dirX[phase_index,ibond] = self.newX[phase_index,j] - self.newX[phase_index,i]
                self.dirY[phase_index,ibond] = self.newY[phase_index,j] - self.newY[phase_index,i]
                self.dirZ[phase_index,ibond] = self.newZ[phase_index,j] - self.newZ[phase_index,i]
            #end for ibond,bond
        # end of phase_index,phase


    def animate(self):
        mlab.clf()
        self.draw_noupdate()
        animator = anim(self.phase_index, self.phase_direction, self.number_of_phase_steps, self.quiverBonds, self.refreshAtoms, self.refreshQuiverBonds,self.refreshBonds)
#        from objbrowser import browse
#        browse(animator)
#        print('In animate called anim')

    def refreshBonds(self, phase_index):
        for ig,glyph in enumerate(self.glyph_bonds):
            glyph.mlab_source.set(x=self.bondx[phase_index,ig], y=self.bondy[phase_index,ig], z=self.bondz[phase_index,ig])
        return

    def refreshQuiverBonds(self, phase_index):
        self.glyph_quiver_bonds.mlab_source.set(x=self.startX[phase_index],
                                                y=self.startY[phase_index],
                                                z=self.startZ[phase_index],
                                                u=self.dirX[phase_index],
                                                v=self.dirY[phase_index],
                                                w=self.dirZ[phase_index])
        return

    def refreshAtoms(self, phase_index):
        for ig,glyph in enumerate(self.glyph_atoms):
            glyph.mlab_source.set(x=self.newX[phase_index,ig], y=self.newY[phase_index,ig], z=self.newZ[phase_index,ig])
        return

    def refreshDisplacements(self):
        mode = self.mode
        self.glyph_displacements.mlab_source.set(u=self.U[mode],v=self.V[mode],w=self.W[mode])
        return

    def draw(self):
        mlab.clf()
        self.draw_noupdate()
        self.drawDisplacements()

    def draw_noupdate(self):
        self.drawUnitCell()
        if self.quiverBonds:
            self.drawQuiverBonds()
        else:
            self.drawBonds()
        self.drawAtoms()
        # self.drawDisplacements()
        return

    def drawDisplacements(self):
        mode = self.mode
        self.glyph_displacements = mlab.quiver3d(self.X,self.Y,self.Z,self.U[mode],self.V[mode],self.W[mode],color=(0.15,0.15,0.15),mode='arrow' )
        return

    def drawUnitCell(self):
        bondx = np.zeros( 2 )
        bondy = np.zeros( 2 )
        bondz = np.zeros( 2 )
        ones  = np.array( (1.0,1.0) )
        for edge in self.cell_edges:
            p1,p2 =edge
            bondx[0] = p1[0]
            bondy[0] = p1[1]
            bondz[0] = p1[2]
            bondx[1] = p2[0]
            bondy[1] = p2[1]
            bondz[1] = p2[2]
            mlab.plot3d(bondx,bondy,bondz,ones, tube_radius=0.15, tube_sides=10, color=(1,0,0) )
        return

    def drawQuiverBonds(self):
        self.glyph_quiver_bonds = mlab.quiver3d(self.startX[0],self.startY[0],self.startZ[0],
                                                self.dirX[0],  self.dirY[0],  self.dirZ[0],
                                                color=(0,0,0), mode='cylinder')

    def drawBonds(self):
        bondx = np.zeros( 2 )
        bondy = np.zeros( 2 )
        bondz = np.zeros( 2 )
        ones  = np.array( (1.0,1.0) )
        self.glyph_bonds = []
        for bond in self.unit_cell.bonds:
            i,j = bond
            bondx[0] = self.X[i]
            bondy[0] = self.Y[i]
            bondz[0] = self.Z[i]
            bondx[1] = self.X[j]
            bondy[1] = self.Y[j]
            bondz[1] = self.Z[j]
            glyph = mlab.plot3d(bondx,bondy,bondz,ones,tube_radius=0.1, tube_sides=10, color=(0,0,0) )
            self.glyph_bonds.append(glyph)
        return

    def drawAtoms(self):
        self.glyph_atoms = []
        for x,y,z,colour,size in zip(self.X, self.Y, self.Z, self.colours, self.radii):
           glyph = mlab.points3d([x],[y],[z],[1],color=colour,resolution=20, scale_factor=size)
           self.glyph_atoms.append(glyph)
        return

    def _get_selected_mode_frequency(self):
        self.mode = self.gui_mode_selection
        self.calculatePhasePositions()
        return self.frequencies[self.gui_mode_selection]

@mlab.animate(delay=20, ui=True, support_movie=False)
def anim(phase_index, phase_direction, number_of_phase_steps, quiverBonds, refreshAtoms, refreshQuiverBonds, refreshBonds ):
    f = mlab.gcf()
    f.scene.anti_aliasing_frames = 0
    while True:
        phase_index = phase_index + phase_direction
        if phase_index >= int(number_of_phase_steps/2):
            phase_direction = -1
        elif phase_index <= -int(number_of_phase_steps/2):
            phase_direction = +1
        f.scene.disable_render = True
        if quiverBonds:
            refreshQuiverBonds(phase_index)
        else:
            refreshBonds(phase_index)
        refreshAtoms(phase_index)
        f.scene.disable_render = False
        yield
    return
