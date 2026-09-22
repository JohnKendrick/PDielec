# The Examples directory for PDielec and PDGui

[Back to PDielec](../README.md)

 | **File**       | **Description**                                         |
 | -------------- | ------------------------------------------------------- |
 | [AbInit](./AbInit/README.md) | AbInit examples  |
 | [ATR](./ATR/README.md) | ATR examples |
 | [Castep](./Castep/README.md) | CASTEP examples  |
 | [Crystal](./Crystal/README.md) | CRYSTAL examples |
 | [Experiment](./Experiment/README.md) | Examples using the Experiment file format |
 | [Gulp](./Gulp/README.md) | GULP examples   |
 | [Helper](./Helper/README.md) | Examples using the API    |
 | [Mie](./Mie/README.md) | Mie scattering examples   |
 | [P2Cif](./P2Cif/README.md) | Generating cif files from the output files    |
 | [Phonopy](./Phonopy/README.md) | PHONOPY examples|
 | [QE](./QE/README.md) | Quantum Espresso examples |
 | [SingleCrystal](./SingleCrystal/README.md) | Examples of single crystal calculations       |
 | [SizeEffects](./SizeEffects/README.md) | Examples including size effects in powder calculations |
 | [Vasp](./Vasp/README.md) | VASP examples   |
 | [VibAnalysis](./VibAnalysis/README.md) | VibAnalysis examples      |

Additional response workflows:

| Directory | Purpose |
| --- | --- |
| [Powder_Raman](./Powder_Raman/README.md) | ZnO powder Raman examples for six input readers |
| [Crystal_Raman](./Crystal_Raman/README.md) | ZnO crystal Raman geometry and coherence examples |
| [Finite_difference](./Finite_difference/README.md) | Infrared and Raman preader recipes using finite-difference JSON |
| [FHI-Aims](./FHI-Aims/README.md) | FHI-Aims response and reader examples |

Run `pdmake command.pdmake` in a directory containing a recipe, or
`pdmake --view command.pdmake` to open a PDGui recipe interactively.
From the repository root, `pdmake tests` runs pytest followed by the configured
example suites. Not every supporting file or notebook is a regression recipe.

## Reading these files locally

Navigation links name `README.md` explicitly. A link to a directory opens a
file manager in many desktop viewers; it does not automatically display the
README as GitHub does. Each example README links back to its parent, and
parent READMEs link to their child directories.

Okular delegates ordinary Markdown links to the desktop's associated
application (see [Okular's link handling](https://github.com/KDE/okular/blob/master/core/document.cpp)).
Opening a README in Okular does not itself make Okular the default Markdown
viewer. If clicking another README opens an editor, select **Okular** as the
default application for Markdown files in your desktop's file associations.
On Linux installations providing `okularApplication_md.desktop`, the equivalent
command is:

```sh
xdg-mime default okularApplication_md.desktop text/markdown
```

Restart Okular after changing the association. KDE and other desktops can use
different association preferences; if the command has no effect, use KDE's
file-association settings for `text/markdown` and move Okular to the top.
Links to data files still open the application associated with that file type.
