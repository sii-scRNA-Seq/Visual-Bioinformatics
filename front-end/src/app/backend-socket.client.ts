import { Inject, Injectable } from '@angular/core';
import { Socket } from 'socket.io-client';

import { Request } from './request';
import { SOCKET } from './socket';

@Injectable({
  providedIn: 'root'
})
export class BackendSocketClient {

  constructor(@Inject(SOCKET) private socket: Socket) {
    console.log('BackendSocketClient connecting to socket');
    this.socket.connect();
    console.log('BackendSocketClient connected to socket');
  }

  sendRequest(request: Request): void {
    this.socket.emit('json', request);
  }

  listen(foo: (res: string) => void): void {
    // if (! this.socket.connected) {
    //   console.log('deferring');
    //   setTimeout(() => this.listen(foo), 500);
    //   return;
    // }
    console.log('BackendSocketClient listening');
    this.socket.on('json', (msg: string) => {
      foo(msg);
    });
  }
}
