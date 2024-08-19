import { InjectionToken, isDevMode } from '@angular/core';
import { Socket, io } from 'socket.io-client';

export const SOCKET = new InjectionToken<Socket>('socket');

const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
export const socket = io(url, {
  autoConnect: false,
  // This is because websocket upgrade always fails in production, and can cause issues with first block loading. May also
  // be fixable by changing nginx configuration
  // See this source: https://socket.io/docs/v4/client-options/#transports and https://socket.io/docs/v4/how-it-works/#upgrade-mechanism
  transports: ['polling'],
});
