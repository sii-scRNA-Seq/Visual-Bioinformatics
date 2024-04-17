import { Injectable } from '@angular/core';
import { BehaviorSubject, Observable, ReplaySubject } from 'rxjs';

import { Block } from './block.interface';
import { Output } from './output';
import { OutputServiceInterface } from './output.service.interface';
import { BackendSocketClientInterface } from './backend-socket.client.interface';

@Injectable({
  providedIn: 'root'
})
export class MockBackendSocketClient implements BackendSocketClientInterface {
  private readonly response$: ReplaySubject<any> = new ReplaySubject<any> (1);
  readonly response: Observable<any> = this.response$.asObservable();
  
  sendRequest(message: any): void {
    const res = {
        text: 'Dummy text'
    };
    this.response$.next(res)
  }
}