/* eslint-disable @typescript-eslint/no-unused-vars */

import { Injectable, isDevMode } from '@angular/core';
import { Socket, io } from 'socket.io-client';

import { BackendSocketClientInterface } from './backend-socket.client.interface';
import { Request } from './request';

@Injectable({
  providedIn: 'root'
})
export class MockBackendSocketClient implements BackendSocketClientInterface {

  private socket: Socket;
  
  constructor() {
    const url = isDevMode() ? 'ws://127.0.0.1:5000' : '';
    // TODO: Are these timeouts required?
    this.socket = io(url, {
      autoConnect: false,
      timeout: 1000000,   // 1000 seconds,
      ackTimeout: 1000000 // 1000 seconds,
    });
    
    this.socket.connect();
  }
  
  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      foo(msg);  
    });
  }

  sendRequest(request: Request): void { }
}