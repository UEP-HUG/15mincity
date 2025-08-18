# Data Input Files - Attribution and Licensing

This folder contains geospatial data files used in the 15-minute city accessibility analysis. All files are used in accordance with their respective licensing terms.

## Administrative Boundaries

### `swissBOUNDARIES3D_1_3_TLM_HOHEITSGEBIET.gpkg`
- **Source**: Swiss Federal Office of Topography (swisstopo)
- **Dataset**: swissBOUNDARIES3D v1.3 - Administrative boundaries (sovereignty areas)
- **License**: Open Government Data License
- **Copyright**: © Swiss Federal Office of Topography (swisstopo)
- **Version**: 1.3
- **Date**: June 11, 2025
- **Attribution Required**: Yes (mandatory source reference)
- **Terms**: https://www.swisstopo.admin.ch/en/terms-of-use-free-geodata-and-geoservices

### `swissBOUNDARIES3D_1_3_TLM_KANTONSGEBIET.gpkg`
- **Source**: Swiss Federal Office of Topography (swisstopo)
- **Dataset**: swissBOUNDARIES3D v1.3 - Cantonal boundaries
- **License**: Open Government Data License  
- **Copyright**: © Swiss Federal Office of Topography (swisstopo)
- **Version**: 1.3
- **Date**: June 11, 2025
- **Attribution Required**: Yes (mandatory source reference)
- **Terms**: https://www.swisstopo.admin.ch/en/terms-of-use-free-geodata-and-geoservices

### `canton_ge.GeoJSON`
- **Source**: Swiss Federal Office of Topography (swisstopo)
- **Dataset**: Derived from swissBOUNDARIES3D - Canton boundaries
- **License**: Open Government Data License
- **Copyright**: © Swiss Federal Office of Topography (swisstopo)
- **Date**: June 2, 2025
- **Attribution Required**: Yes (per swisstopo terms of use)
- **Terms**: https://www.swisstopo.admin.ch/en/terms-of-use-free-geodata-and-geoservices

### `communes_ge.GeoJSON`
- **Source**: Swiss Federal Office of Topography (swisstopo)  
- **Dataset**: Derived from swissBOUNDARIES3D - Municipality boundaries
- **License**: Open Government Data License
- **Copyright**: © Swiss Federal Office of Topography (swisstopo)
- **Date**: Current (updated regularly)
- **Attribution Required**: Yes (per swisstopo terms of use)
- **Terms**: https://www.swisstopo.admin.ch/en/terms-of-use-free-geodata-and-geoservices

## Water Bodies

### `GEO_LAC.gpkg`
- **Source**: SITG (Système d'information du territoire à Genève)
- **Dataset**: Hydrographic features - Lake Geneva, Rhône and Arve rivers
- **Description**: Surface hydrographic layer for cartographic production representing Lake Geneva, Rhône and Arve rivers
- **Based on**: Survey data from Geneva, Vaud and Valais cantons + French IGN data
- **License**: Open Government Data License
- **Copyright**: © SITG
- **Update Frequency**: Regular updates
- **Source URL**: https://sitg.ge.ch/donnees/geo-lac
- **Attribution Required**: Yes

## Population Data

### `GML_AGGLO_CARREAU_200/` (folder)
- **Source**: SITG (Système d'information du territoire à Genève)
- **Dataset**: Greater Geneva Population Grid (200m resolution)
- **License**: Open Government Data License
- **Data Vintage**: 2019
- **Date Retrieved**: October 14, 2024
- **Coverage**: Cross-border Greater Geneva metropolitan area
- **Underlying Sources**: 
  - Switzerland: STATPOP (© Swiss Federal Statistical Office, 2019)
  - France: FiLoSoFi (© INSEE, 2019)
- **Access**: https://sitg.ge.ch/donnees/agglo-carreau-population
- **Attribution**: SITG, SFSO, INSEE

## Data Processing Notes

1. **File Formats**: Original data converted to standardized formats (GeoJSON/GPKG) for analysis
2. **Coordinate System**: All data projected to appropriate coordinate reference system for analysis
3. **Cross-Border Integration**: Administrative boundaries and population data harmonized across Swiss-French border
4. **Attribution Compliance**: All files used comply with respective licensing terms requiring source attribution

## License Summary

All geospatial data files are available under **Open Government Data licenses** requiring:
- Source attribution in publications
- Copyright notice preservation  
- License terms acknowledgment
- No commercial restrictions
- No modification restrictions

## Citation Requirements

When using these data files, include appropriate attribution:
- **swisstopo data**: "© Swiss Federal Office of Topography (swisstopo), [year]"
- **SITG data**: "© SITG, [year]"
- **Combined datasets**: Include all relevant copyright holders

For detailed citation requirements, see main project README.md and original data source documentation.