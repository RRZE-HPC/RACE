/*
 * =======================================================================================
 *
 *   RACE: Recursicve Algebraic Coloring Engine
 *   Copyright (C) 2019, RRZE, Friedrich-Alexander-Universität Erlangen-Nürnberg
 *   Author: Christie Alappat
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Affero General Public License as
 *   published by the Free Software Foundation, either version 3 of the
 *   License, or (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Affero General Public License for more details.
 *
 *   You should have received a copy of the GNU Affero General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 * =======================================================================================
 */

#include "signal.h"

Signal::Signal(int RACE_BLOCKCTR)
{
    lock = new pthread_mutex_t;
    pthread_mutex_init(lock, NULL);
    signal = new pthread_cond_t;
    pthread_cond_init(signal, NULL);
    preSignal = new spin_cond_t;
    spin_cond_init(preSignal, RACE_BLOCKCTR);
}

Signal::~Signal()
{
    delete lock;
    delete signal;
    delete preSignal;
}
