Version 10 migration
====================

Version 10 adds powder and crystal Raman scenarios. Before reusing an older
script or saved scenario, check the following changes.

Python scripts
--------------

Public methods and functions formerly written in camelCase have been renamed
to snake_case throughout the numerical modules, readers and GUI. Older scripts
calling those names need updating. Use the current API reference and the
``script.py`` files in the examples as templates; do not assume every rename
can be obtained by a mechanical text replacement. Save a new settings script
from the GUI after checking an older calculation.

Crystal Raman geometry
----------------------

``Collection angle = -1`` selects automatic collection. On the superstrate
side this now means retro-backscattering; on the substrate side it means
collinear forward scattering. Saved automatic scenarios can therefore change
at nonzero incidence. Normal-incidence automatic scenarios are unchanged.
To retain specular reflection, explicitly set the collection angle equal to
the incidence angle. Other negative collection angles are explicit signed
angles, not automatic-selection sentinels.

``Approximate ES`` defaults to False. When enabled, it evaluates the reciprocal
detector field once at the laser frequency instead of at each Stokes frequency.
It reuses the incident field only when the detector optical system and signed
angle match. Review these settings before comparing old and new spectra;
see :doc:`CrystalRaman` and :doc:`settings_reference`.

Spreadsheets and regression references
--------------------------------------

Settings-sheet Raman activities now follow the activity units selected in the
Settings frequency table, with units in the column headers. ``pdmake`` checks
this sheet as well as result sheets and reports a checked sheet missing from
one workbook. Review differences before regenerating reference workbooks.
See :doc:`software` for the conditional Raman and infrared output sheets.
