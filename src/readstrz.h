/**
 * @file readstrz.h
 *
 * @brief Conversion of files in Strahlenzentrum (STRZ) data format
 *
 * @author Stefan Schippers
 * @verbatim
// SPDX-License-Identifier: MIT
   @endverbatim
 *
 */
#pragma once

namespace STRZ {

/////////////////////////////////////////////////////////////////////////////////
/**
 * @brief Reads binary data files from Giessen electron-ion crossed-beams setup
 * and converts to ascii
 */
int ConvertStrzDataFile(void);

} // end namespace STRZ
