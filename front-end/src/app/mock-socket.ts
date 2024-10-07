/* eslint-disable @typescript-eslint/no-unused-vars */

import { Injectable } from '@angular/core';
import { Request } from './request';

@Injectable({
  providedIn: 'root'
})
export class MockSocket {

  connect(): void { }

  emit(eventName: string, request: Request): void { }

  on(eventName: string, callback: (msg: string) => void): void { }
}
