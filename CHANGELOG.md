# Changelog

## [0.1.4] - 2020-09-25

### Added
- *pyBLOCHUS* - core functionality of *BLOCHUS* as python module without a graphical user interface.

### Changed
- Minor GUI improvements and consistency clean ups.

### Fixed
- bug inside the pulse modulation function

## [0.1.3] - 2020-05-22

### Added
- Two additional examples that demonstrate the creation of lookup tables (one example for pre-polarization switch-off ramps and one example for adiabatic excitation pulses).

### Changed
- Minor GUI improvements and consistency clean ups.

## [0.1.2] - 2020-05-17

### Changed
- Updated `README.md` and `CHANGELOG.md`.

### Fixed

- When switching the *Pulse* panel on/off, the quality factor settings did not behave as expected.

## [0.1.1] - 2020-05-16

### Added
- An animation of the latest result can be played for all parameter combinations (e.g. *Pulse + Relaxation* or *Pre-polarization switch-off + Pulse + Relaxation*, etc.).

### Fixed

- The adiabatic quality p was incorrectly calculated in `onPushRun.m` when the total simulation time exceeded the actual switch-off ramp time.
- When exporting the *Current View* figure from the menu, both axes for the pulse setup (frequency and current modulation) are now correctly shown.
- Fixed typos in comments and `README.md`.

## [0.1.0] - 2020-05-15

Initial Version

[0.1.4]: https://github.com/ThoHiller/nmr-blochus/compare/v0.1.3...v0.1.4
[0.1.3]: https://github.com/ThoHiller/nmr-blochus/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/ThoHiller/nmr-blochus/compare/v.0.1.1...v0.1.2
[0.1.1]: https://github.com/ThoHiller/nmr-blochus/compare/v.0.1.0...v.0.1.1
[0.1.0]: https://github.com/ThoHiller/nmr-blochus/releases/tag/v.0.1.0
