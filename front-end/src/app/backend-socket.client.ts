import { Inject, Injectable } from '@angular/core';
import { Socket } from 'socket.io-client';

import { Request } from './request';
import { SOCKET } from './socket';

@Injectable({
  providedIn: 'root'
})
export class BackendSocketClient {

  constructor(@Inject(SOCKET) private socket: Socket) {
    this.socket.connect();
  }

  /**
   * Sends a WebSocket request to the backend.
   *
   * @param request - The request which will be sent to the backend
   */
  sendRequest(request: Request): void {
    this.socket.emit('json', request);
  }

  /**
   * Listens for WebSocket responses from the backend and executes a function on any responses.
   *
   * @param foo - The function to be executed on each response
   */
  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      foo(msg);
    });
  }
}
