/**
 * @file readstrz.h
 *
 * @brief Conversion of files in Strahlenzentrum (STRZ) data format
 *
 * @author Stefan Schippers
 * @verbatim
   $Id: readstrz.h 2037 2026-07-17 15:21:56Z iamp $
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
