import { InjectionToken, isDevMode } from '@angular/core';
import { Socket, io } from 'socket.io-client';

export const SOCKET = new InjectionToken<Socket>('socket');
const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
// TODO: Are these timeouts required?
export const socket = io(url, {
  autoConnect: false,
  timeout: 1000000,   // 1000 seconds,
  ackTimeout: 1000000 // 1000 seconds,
});