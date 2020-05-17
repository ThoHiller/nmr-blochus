# Changelog

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

[0.1.2]: https://github.com/ThoHiller/nmr-blochus/compare/v.0.1.1...v0.1.2
[0.1.1]: https://github.com/ThoHiller/nmr-blochus/compare/v.0.1.0...v.0.1.1
[0.1.0]: https://github.com/ThoHiller/nmr-blochus/releases/tag/v.0.1.0
