import { InjectionToken, isDevMode } from '@angular/core';
import { Socket, io } from 'socket.io-client';

export const SOCKET = new InjectionToken<Socket>('socket');

const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
export const socket = io(url, { autoConnect: false });
