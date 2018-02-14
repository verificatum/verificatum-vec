
/*
 * Copyright 2008-2018 Douglas Wikstrom
 *
 * This file is part of Verificatum Elliptic Curve library (VEC).
 *
 * VEC is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * VEC is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General
 * Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with VEC. If not, see
 * <http://www.gnu.org/licenses/>.
 */

#ifndef TEMPLATES_H
#define TEMPLATES_H

#define CAT(X,Y) X##Y
#define FUNCTION_NAME(X,Y) CAT(X,Y)

/*
 * We use compiler flags that enforce that unused variables are
 * flagged as errors. We need this in some cases where we are forced
 * to keep parameters to satisfy a given parameter signature.
 */
#define VEC_UNUSED(x) ((void)(x))

#endif /* TEMPLATES_H */
