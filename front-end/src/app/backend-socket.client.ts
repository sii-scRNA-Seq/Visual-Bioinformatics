import { Injectable, isDevMode } from '@angular/core';
import { Observable, ReplaySubject } from 'rxjs';
import { Socket, io } from 'socket.io-client';

import { BackendSocketClientInterface } from './backend-socket.client.interface';

@Injectable({
  providedIn: 'root'
})
export class BackendSocketClient implements BackendSocketClientInterface {

  private socket: Socket;

  private readonly response$: ReplaySubject<any> = new ReplaySubject<any> (1);
  readonly response: Observable<any> = this.response$.asObservable();

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

  sendRequest(message: any): void {
    this.socket.emit('json', message);
  }

  listen(foo: (res: string) => void): void {
    this.socket.on('json', (msg: string) => {
      foo(msg);  
    });
  }
}